"""Main orchestrator (CLAUDE.md Step 8).

Fetches PubMed / UniProt / OpenTargets for each gene, extracts structured
annotations, merges them, enriches with STRING / GTEx / CellxGene, then builds
the target network and prioritized table.

The per-gene work is split into composable stages — fetch → extract → merge →
enrich — each a small async function taking ``(gene, …, config)``. ``run_gene``
chains them; ``main`` runs ``run_gene`` across all genes with a concurrency cap
and per-gene error isolation.

Run:  python pipeline.py
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import logging
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter
from typing import Any

import httpx
from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

from src.config import config, PipelineConfig
from src.utils import setup_logging, load_gene_list, load_ref_set
from src.fetchers.pubmed import search_pmids, fetch_abstracts
from src.fetchers.uniprot import fetch_uniprot
from src.fetchers.opentargets import fetch_opentargets
from src.fetchers.string_db import fetch_string
from src.fetchers.cellxgene import fetch_cellxgene
from src.extractor import (
    extract_from_text,
    extract_from_uniprot,
    extract_from_opentargets,
)
from src.merger import merge_annotations
from src.filters.gtex_safety import assess_safety
from src.network import (
    build_target_network,
    compute_priority_scores,
    save_network,
    save_prioritized_tsv,
)

ANNOTATIONS_JSONL = Path("outputs/annotations.jsonl")

log = logging.getLogger("bio_annot.pipeline")

# Canonical Reactome reference, loaded once and cached. The merge stage needs it,
# but the stage signatures are (gene, …, config), so it is sourced here rather
# than threaded through every call.
_REACTOME_REF_CACHE: set[str] | None = None


def _reactome_ref() -> set[str]:
    """Lazily load and cache the canonical Reactome name set."""
    global _REACTOME_REF_CACHE
    if _REACTOME_REF_CACHE is None:
        _REACTOME_REF_CACHE = load_ref_set("refs/reactome_pathways.txt")
    return _REACTOME_REF_CACHE


# Run-report box rendering. _BOX_WIDTH is the inner width between the ║ borders;
# _row pads/truncates content to exactly that width so every right border lines
# up, and _divider draws a horizontal rule with the given corner/junction chars.
_BOX_WIDTH = 54


def _divider(left: str = "╠", right: str = "╣") -> str:
    """A horizontal rule: ``left`` + box-width of ═ + ``right``."""
    return f"{left}{'═' * _BOX_WIDTH}{right}"


def _row(text: str = "") -> str:
    """A box row: content padded/truncated to the box width, framed by ║."""
    return f"║{text[:_BOX_WIDTH]:<{_BOX_WIDTH}}║"


@dataclass
class RunStats:
    """Aggregate run metrics: gene outcomes, LLM token usage, cost, runtime."""
    start_time: float = field(default_factory=perf_counter)
    genes_total: int = 0
    genes_succeeded: int = 0
    genes_failed: int = 0
    genes_cached: int = 0
    llm_calls: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    cache_read_tokens: int = 0
    cache_created_tokens: int = 0
    noncanonical_total: int = 0
    noncanonical_flagged: list[str] = field(default_factory=list)

    def add_usage(self, usage: dict) -> None:
        self.llm_calls += 1
        self.input_tokens += usage.get("input_tokens", 0)
        self.output_tokens += usage.get("output_tokens", 0)
        self.cache_read_tokens += usage.get(
            "cache_read_input_tokens", 0)
        self.cache_created_tokens += usage.get(
            "cache_creation_input_tokens", 0)

    @property
    def runtime_seconds(self) -> float:
        return perf_counter() - self.start_time

    @property
    def cache_hit_rate(self) -> float:
        if self.llm_calls == 0:
            return 0.0
        return self.cache_read_tokens / max(
            self.cache_read_tokens + self.cache_created_tokens, 1)

    @property
    def estimated_cost_usd(self) -> float:
        # Sonnet 4.6 pricing: $3/$15 per 1M in/out
        # Cache read: $0.30/1M, cache write: $3.75/1M
        input_cost = (self.input_tokens / 1_000_000) * 3.0
        output_cost = (self.output_tokens / 1_000_000) * 15.0
        cache_read_cost = (self.cache_read_tokens / 1_000_000) * 0.30
        cache_write_cost = (
            self.cache_created_tokens / 1_000_000) * 3.75
        return input_cost + output_cost + cache_read_cost + cache_write_cost

    def print_report(self, annotations: dict) -> None:
        # Count NON-CANONICAL across all annotations
        for gene, ann in annotations.items():
            for p in ann.get("pathways", []):
                if p.startswith("NON-CANONICAL:"):
                    self.noncanonical_total += 1
                    self.noncanonical_flagged.append(
                        f"{gene}: {p[14:].strip()}")

        total_pathways = sum(
            len(ann.get("pathways", []))
            for ann in annotations.values()
        )
        nc_pct = self.noncanonical_total / max(total_pathways, 1) * 100

        mins = int(self.runtime_seconds // 60)
        secs = int(self.runtime_seconds % 60)

        lines = [
            "",
            _divider("╔", "╗"),
            _row(f"{'BioAnnot Run Report':^{_BOX_WIDTH}}"),
            _divider(),
            _row(" Genes"),
            _row(f"   Total:      {self.genes_total:<4}  Succeeded: {self.genes_succeeded:<4}"),
            _row(f"   Failed:     {self.genes_failed:<4}  Cached:    {self.genes_cached:<4}"),
            _divider(),
            _row(" Pathway Quality"),
            _row(f"   NON-CANONICAL: {self.noncanonical_total}/{total_pathways} ({nc_pct:.1f}%)"),
            _divider(),
            _row(" LLM Usage"),
            _row(f"   Calls:      {self.llm_calls}"),
            _row(f"   Input:      {self.input_tokens:<8} tokens"),
            _row(f"   Output:     {self.output_tokens:<8} tokens"),
            _row(f"   Cache read: {self.cache_read_tokens:<8} tokens"),
            _row(f"   Cache hit:  {self.cache_hit_rate*100:.1f}%"),
            _row(f"   Est. cost:  ${self.estimated_cost_usd:.4f}"),
            _divider(),
            _row(f" Runtime: {mins}m {secs}s"),
            _divider("╚", "╝"),
        ]
        report = "\n".join(lines)
        log.info(report)
        print(report)


async def run_fetch_stage(gene: str, config: PipelineConfig) -> dict[str, Any]:
    """Fetch raw data from PubMed, UniProt, OpenTargets.

    Returns a dict with the assembled PubMed abstract text (plus its PMIDs) and
    the raw UniProt / OpenTargets records (``{}`` when a source has no entry).
    """
    pmids = await search_pmids(gene, config.pubmed_max_results)
    abstracts = await fetch_abstracts(pmids)
    pubmed_text = "\n\n".join(
        f"PMID:{a['pmid']}\n{a['abstract']}" for a in abstracts
    )
    up_data = await fetch_uniprot(gene)
    ot_data = await fetch_opentargets(gene)
    return {
        "pmids": pmids,
        "pubmed_text": pubmed_text,
        "uniprot": up_data,
        "opentargets": ot_data,
    }


async def run_extract_stage(
    gene: str, fetched: dict[str, Any], config: PipelineConfig
) -> tuple[list[dict], list[dict]]:
    """Extract structured annotations from each source.

    Runs the LLM extractor over PubMed text, UniProt, and OpenTargets, persists
    the (unfiltered) per-source extractions to ``outputs/raw/{gene}_raw.json``,
    and returns ``(sources, usages)`` — the sources clearing the confidence gate,
    plus the per-call token-usage dicts for stats accounting.
    """
    usages: list[dict] = []

    pubmed_ann, pubmed_usage = await extract_from_text(
        gene, fetched["pubmed_text"], fetched["pmids"]
    )
    usages.append(pubmed_usage)

    up_data = fetched["uniprot"]
    ot_data = fetched["opentargets"]

    uniprot_ann = None
    if up_data:
        uniprot_ann, uniprot_usage = await extract_from_uniprot(gene, up_data)
        usages.append(uniprot_usage)

    ot_ann = None
    if ot_data:
        ot_ann, ot_usage = await extract_from_opentargets(gene, ot_data)
        usages.append(ot_usage)

    # Persist the raw per-source LLM extractions for provenance/debugging.
    # (STRING partners are persisted with the merged record, not here.)
    raw_path = Path("outputs/raw") / f"{gene}_raw.json"
    with open(raw_path, "w", encoding="utf-8") as f:
        json.dump(
            {"pubmed": pubmed_ann, "uniprot": uniprot_ann, "opentargets": ot_ann},
            f,
            indent=2,
        )

    # Keep only high-confidence sources for the merge.
    sources = [
        a
        for a in [pubmed_ann, uniprot_ann, ot_ann]
        if a and a.get("confidence", 0) >= config.confidence_threshold
    ]
    return sources, usages


async def run_merge_stage(
    gene: str, extractions: list[dict], config: PipelineConfig
) -> tuple[dict | None, dict | None]:
    """Merge and resolve conflicts across sources.

    Returns ``(merged, usage)``. ``merged`` is ``None`` (a clean skip, logged)
    when no source cleared the confidence gate, so the gene is excluded from
    output without being treated as a failure. ``usage`` is the merge call's
    token-usage dict, or ``None`` when no LLM merge call was made (single source
    or no sources).
    """
    if not extractions:
        log.warning("No high-confidence sources for %s", gene)
        return None, None
    return await merge_annotations(gene, extractions, _reactome_ref())


async def run_enrich_stage(
    gene: str, merged: dict, config: PipelineConfig
) -> dict:
    """Enrich the merged record with STRING, GTEx, and CellxGene.

    Mutates and returns ``merged`` with ``string_interactors``,
    ``safety_assessment``, and (when enabled) ``cellxgene_expression`` plus
    measured cell types unioned into ``cellular_states``.
    """
    # STRING PPI partners (factual, not LLM-extracted) for the network builder.
    string_interactors = await fetch_string(gene)
    merged["string_interactors"] = string_interactors

    # GTEx normal-tissue safety assessment for the network scorer.
    safety = assess_safety(gene)
    merged["safety_assessment"] = safety
    if safety.get("safety_flag"):
        log.warning(
            "Safety flag for %s: high normal expression in %d sensitive tissues %s "
            "(max %.1f TPM)",
            gene,
            safety.get("tissue_count_above_threshold", 0),
            safety.get("high_expression_tissues", []),
            safety.get("max_tpm", 0.0),
        )

    # CellxGene Census single-cell expression — grounds cellular_states in
    # measured per-cell-type expression for the configured tissue.
    if config.enable_cellxgene:
        census_data = await fetch_cellxgene(gene)
        merged["cellxgene_expression"] = {
            "tissue": config.census_tissue,
            "census_version": config.census_version,
            "top_cell_types": [
                {"cell_type": k, "mean_expr": v}
                for k, v in list(census_data.items())[:10]
            ],
            "cell_type_count": len(census_data),
        }
        # Union the top 5 measured cell types into cellular_states (prefixed so
        # their provenance is distinguishable from LLM-extracted states), order-
        # preserving and deduped.
        top5 = [f"CellxGene: {ct}" for ct in list(census_data.keys())[:5]]
        merged["cellular_states"] = list(
            dict.fromkeys(merged.get("cellular_states", []) + top5)
        )

    return merged


def make_cache_key(gene: str, config: PipelineConfig) -> str:
    """Content-addressed cache key for a gene's result.

    The key digests the inputs that determine the output — gene, disease
    context, the extraction/merge models, the PubMed depth, and the CellxGene
    tissue/version — so changing any of them invalidates the cache for that gene.
    """
    components = "|".join([
        gene,
        config.disease_context,
        config.extraction_model,
        config.merge_model,
        str(config.pubmed_max_results),
        config.census_tissue,
        config.census_version,
    ])
    return hashlib.md5(components.encode()).hexdigest()[:12]


def read_cache(gene: str, config: PipelineConfig) -> dict | None:
    """Return a gene's cached result, or None on miss / cache disabled / force.

    Honors ENABLE_CACHE (off → no cache) and FORCE_RERUN (on → ignore existing
    cache so all stages recompute; the fresh result is still written back).
    """
    if not config.enable_cache or config.force_rerun:
        return None
    path = Path(config.cache_dir) / f"{gene}_{make_cache_key(gene, config)}.json"
    if path.exists():
        log.info("Gene %s: cache hit → skipping all stages", gene)
        return json.loads(path.read_text())
    return None


def write_cache(gene: str, config: PipelineConfig, result: dict) -> None:
    """Persist a gene's result to the resume cache (no-op if caching disabled)."""
    if not config.enable_cache:
        return
    Path(config.cache_dir).mkdir(parents=True, exist_ok=True)
    path = Path(config.cache_dir) / f"{gene}_{make_cache_key(gene, config)}.json"
    path.write_text(json.dumps(result, indent=2, default=str))
    log.info("Gene %s: cached to %s", gene, path)


async def run_gene(
    gene: str, config: PipelineConfig, stats: RunStats
) -> dict | None:
    """Run all stages for a single gene, accumulating LLM usage into ``stats``.

    Returns the enriched annotation, or ``None`` when no high-confidence source
    survived the merge gate. On a resume-cache hit, returns the cached result
    without running any stage (and counts it in ``stats.genes_cached``).
    """
    # Check the on-disk resume cache first.
    cached = read_cache(gene, config)
    if cached is not None:
        stats.genes_cached += 1
        # Marker so the progress bar can show "cached"; stripped in main() before
        # the record is written to final_annotations.json.
        cached["_from_cache"] = True
        return cached

    fetched = await run_fetch_stage(gene, config)
    extractions, extract_usages = await run_extract_stage(gene, fetched, config)
    for usage in extract_usages:
        stats.add_usage(usage)

    merged, merge_usage = await run_merge_stage(gene, extractions, config)
    if merged is None:
        return None
    if merge_usage:
        stats.add_usage(merge_usage)

    enriched = await run_enrich_stage(gene, merged, config)

    # Append the merged record to the run's annotations log, then cache it.
    with open(ANNOTATIONS_JSONL, "a", encoding="utf-8") as f:
        f.write(json.dumps(enriched) + "\n")
    write_cache(gene, config, enriched)
    return enriched


# Transient transport-level errors worth retrying. Logic/data errors
# (ValueError, KeyError, …) are intentionally excluded so they fail fast.
RETRYABLE_ERRORS = (
    httpx.ConnectTimeout,
    httpx.ConnectError,
    httpx.ReadTimeout,
    asyncio.TimeoutError,
)


async def run_gene_with_retry(
    gene: str, config: PipelineConfig, stats: RunStats
) -> dict | None:
    """Run all stages for a gene with retry on transient errors.

    Up to 3 attempts with exponential backoff (2s, 4s) on transport errors;
    non-retryable errors propagate immediately.
    """
    for attempt in range(1, 4):  # 3 attempts
        try:
            return await run_gene(gene, config, stats)
        except RETRYABLE_ERRORS as e:
            if attempt == 3:
                log.error(
                    "Gene %s failed after 3 attempts: %s",
                    gene, e
                )
                raise
            wait = 2 ** attempt  # 2s, 4s
            log.warning(
                "Gene %s attempt %d failed (%s), retrying in %ds",
                gene, attempt, type(e).__name__, wait
            )
            await asyncio.sleep(wait)


async def main() -> None:
    setup_logging(config.log_level)
    Path("outputs/raw").mkdir(parents=True, exist_ok=True)
    # Start each run with a fresh annotations log so re-runs don't accumulate.
    ANNOTATIONS_JSONL.unlink(missing_ok=True)

    log.info("Running in disease context: %s", config.disease_context)
    genes = load_gene_list("inputs/target_genes.txt")
    log.info("Processing %d genes", len(genes))

    stats = RunStats()
    stats.genes_total = len(genes)

    # Process genes with a concurrency limit (respect API rate limits).
    sem = asyncio.Semaphore(config.semaphore_limit)

    pbar = tqdm(
        total=len(genes),
        desc="BioAnnot",
        unit="gene",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} "
                   "[{elapsed}<{remaining}, {rate_fmt}]"
    )

    async def bounded_with_progress(gene: str) -> dict | None:
        async with sem:
            pbar.set_postfix(gene=gene, stage="running", refresh=True)
            result = await run_gene_with_retry(gene, config, stats)
            status = "cached" if isinstance(result, dict) and \
                     result.get("_from_cache") else "done"
            pbar.set_postfix(gene=gene, stage=status, refresh=True)
            pbar.update(1)
            return result

    tasks = [bounded_with_progress(g) for g in genes]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    pbar.close()

    # Error isolation: one gene's failure (e.g. a ConnectTimeout) must not kill
    # the whole run — log it, exclude it, and keep the rest.
    final_annotations: dict = {}
    failed = []
    for gene, result in zip(genes, results):
        if isinstance(result, Exception):
            log.error("Gene %s failed: %s, skipping", gene, result)
            failed.append(gene)
            stats.genes_failed += 1
        elif result:
            # Strip the progress-bar cache marker so it never reaches output.
            result.pop("_from_cache", None)
            final_annotations[gene] = result
            stats.genes_succeeded += 1
    if failed:
        log.warning("Failed genes (excluded from output): %s", failed)

    # Write final merged JSON
    with open("outputs/final_annotations.json", "w", encoding="utf-8") as f:
        json.dump(final_annotations, f, indent=2)
    log.info(
        "Wrote %d annotations → outputs/final_annotations.json", len(final_annotations)
    )

    # Build network and prioritize
    G = build_target_network(final_annotations)
    save_network(G, "outputs/target_network.gpickle")
    scores = compute_priority_scores(G, config.disease_context)
    save_prioritized_tsv(scores, "outputs/prioritized_targets.tsv")
    log.info("Top 5 targets: %s", [s["gene"] for s in scores[:5]])

    # Optionally refresh the pathway synonym map from this run's NON-CANONICAL
    # names, so the next run's fuzzy canonicalization picks them up.
    if config.auto_update_synonyms:
        log.info("AUTO_UPDATE_SYNONYMS=true — updating refs/pathway_synonyms.json")
        subprocess.run([sys.executable, "scripts/build_synonyms.py"], check=False)

    # Run report: gene outcomes, pathway quality, LLM usage, cost, runtime.
    stats.print_report(final_annotations)


if __name__ == "__main__":
    asyncio.run(main())
