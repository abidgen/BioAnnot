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
from src.pathways import SYNONYMS_PATH
from src.filters.gtex_safety import assess_safety
from src.network import (
    build_target_network,
    compute_priority_scores,
    save_network,
    save_prioritized_tsv,
)

ANNOTATIONS_JSONL = Path("outputs/annotations.jsonl")

# Canonical Reactome reference path — used both to load the name set and to
# fingerprint the file into the final-layer cache key.
REACTOME_REF_PATH = "refs/reactome_pathways.txt"

log = logging.getLogger("bio_annot.pipeline")

# Canonical Reactome reference, loaded once and cached. The merge stage needs it,
# but the stage signatures are (gene, …, config), so it is sourced here rather
# than threaded through every call.
_REACTOME_REF_CACHE: set[str] | None = None


def _reactome_ref() -> set[str]:
    """Lazily load and cache the canonical Reactome name set."""
    global _REACTOME_REF_CACHE
    if _REACTOME_REF_CACHE is None:
        _REACTOME_REF_CACHE = load_ref_set(REACTOME_REF_PATH)
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
    genes_remerged: int = 0
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
            _row(f"   Remerged:   {self.genes_remerged:<4} (raw-cache, no fetch/extract)"),
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
    the per-source extractions to ``outputs/raw/{gene}_raw.json``, and returns
    ``(extractions, usages)`` — every successful per-source extraction
    (**pre-confidence-filter**) plus the per-call token-usage dicts for stats
    accounting. The confidence gate is applied later, at consumption in
    ``run_gene``, so the raw cache stores the unfiltered extractions and a
    ``CONFIDENCE_THRESHOLD`` change is honored by re-filtering on a re-merge
    rather than forcing re-extraction.
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

    # Return every successful per-source extraction, unfiltered; run_gene applies
    # the confidence gate at consumption (see docstring).
    extractions = [a for a in [pubmed_ann, uniprot_ann, ot_ann] if a]
    return extractions, usages


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


# ── Two-layer resume cache ────────────────────────────────────────────────────
#
# The cache is split into two independent layers so that editing the pathway
# synonym map or the Reactome reference can refresh annotations WITHOUT paying for
# fetch + extract again:
#
#   Layer 1 — raw extraction cache   outputs/cache/raw/{gene}_{extract_key}.json
#       The unfiltered per-source extractions (the confidence gate is applied at
#       merge time, not here). Keyed by fetch + extract inputs ONLY (gene, disease
#       context, extraction model, PubMed depth, extraction-prompt version).
#       Excludes the synonym map, the Reactome reference, the confidence threshold,
#       and all merge/enrich params — so editing any of those never invalidates it.
#
#   Layer 2 — final enriched cache   outputs/cache/final/{gene}_{full_key}.json
#       The final enriched record (unchanged contents/schema). Keyed by the
#       extract_key PLUS the synonym-map and Reactome-reference file hashes and the
#       downstream config (confidence threshold, merge model, merge/enrich params)
#       — so editing any of those invalidates only this layer, and the raw layer
#       (hence fetch + extract) stays valid for a fetch/extract-free re-merge.
#
# Bump EXTRACT_PROMPT_VERSION whenever the extractor's prompt or output contract
# changes in a way that should invalidate cached raw extractions.
EXTRACT_PROMPT_VERSION = "1"

RAW_CACHE_SUBDIR = "raw"
FINAL_CACHE_SUBDIR = "final"


# Memoize file fingerprints within a process: (path, mtime_ns, size) → digest.
# The synonym map and Reactome reference don't change mid-run, so each is read +
# hashed at most once even though make_full_key is called ~2x per gene. Keying on
# stat() means an edit (new mtime/size) still produces a fresh digest.
_FINGERPRINT_CACHE: dict[tuple, str] = {}


def _file_fingerprint(path) -> str:
    """A short content hash of a file, or "absent" if it does not exist.

    Used to fold the synonym map and Reactome reference into the final cache key
    so any edit to either file changes the key (and only the final layer's key).
    The content hash is computed once per (path, mtime, size) and cached, so the
    large Reactome file is not re-read on every key computation in a run.
    """
    p = Path(path)
    try:
        st = p.stat()
    except OSError:
        return "absent"
    cache_key = (str(p), st.st_mtime_ns, st.st_size)
    digest = _FINGERPRINT_CACHE.get(cache_key)
    if digest is None:
        digest = hashlib.md5(p.read_bytes()).hexdigest()[:12]
        _FINGERPRINT_CACHE[cache_key] = digest
    return digest


def make_extract_key(gene: str, config: PipelineConfig) -> str:
    """Cache key for the raw extraction layer (fetch + extract inputs only).

    Digests the gene, the disease context (which shapes the PubMed query and the
    extractor system prompt), the extraction model, the PubMed depth, and the
    extraction-prompt version. Deliberately excludes the synonym map, the Reactome
    reference, and every merge/enrich parameter so synonym/reference edits leave
    the raw layer valid.

    DO NOT add the synonym-map or Reactome-reference hash here — those belong in
    ``make_full_key`` (the final layer). Folding them in would couple the two cache
    layers and silently break the fetch/extract-free re-merge after a synonym
    edit. The test
    ``tests/test_cache.py::test_extract_key_excludes_merge_and_enrich_params``
    guards this invariant.
    """
    components = "|".join([
        gene,
        config.disease_context,
        config.extraction_model,
        str(config.pubmed_max_results),
        EXTRACT_PROMPT_VERSION,
    ])
    return hashlib.md5(components.encode()).hexdigest()[:12]


# Config fields that change the FINAL record given the same raw extractions, i.e.
# everything consumed downstream of extraction: the confidence gate (applied at
# merge time), the merge model, and every merge/enrich knob. Keep this in sync
# with run_gene's confidence filter + run_merge_stage + run_enrich_stage — adding
# a merge/enrich parameter means adding it here so the final cache invalidates
# when it changes. (The synonym-map and Reactome-reference file hashes are folded
# in separately, since they are file content rather than config fields.)
_FINAL_KEY_CONFIG_FIELDS = (
    "confidence_threshold",
    "merge_model",
    "census_tissue",
    "census_version",
    "enable_cellxgene",
    "string_min_score",
    "string_limit",
    "gtex_tpm_threshold",
    "gtex_min_tissues",
)


def make_full_key(gene: str, config: PipelineConfig) -> str:
    """Cache key for the final enriched layer (the whole chain).

    Builds on ``make_extract_key`` and adds everything that can change the final
    record given the same raw extractions: the synonym-map and Reactome-reference
    file fingerprints plus every field in ``_FINAL_KEY_CONFIG_FIELDS`` (confidence
    gate, merge model, and the merge/enrich params). Editing the synonym map or the
    Reactome reference changes only this key, never the extract key.
    """
    components = [
        make_extract_key(gene, config),
        _file_fingerprint(SYNONYMS_PATH),
        _file_fingerprint(REACTOME_REF_PATH),
    ]
    components += [str(getattr(config, field)) for field in _FINAL_KEY_CONFIG_FIELDS]
    return hashlib.md5("|".join(components).encode()).hexdigest()[:12]


# Backward-compatible alias: the historical single-cache key is the final-layer
# key. Existing callers and tests that import ``make_cache_key`` keep working.
make_cache_key = make_full_key


def _raw_cache_path(gene: str, config: PipelineConfig) -> Path:
    return (Path(config.cache_dir) / RAW_CACHE_SUBDIR
            / f"{gene}_{make_extract_key(gene, config)}.json")


def _final_cache_path(gene: str, config: PipelineConfig) -> Path:
    return (Path(config.cache_dir) / FINAL_CACHE_SUBDIR
            / f"{gene}_{make_full_key(gene, config)}.json")


def _load_cache_file(path: Path) -> dict | None:
    """Load a JSON cache file; treat a corrupt/unreadable one as a miss.

    A truncated file (a run interrupted mid-write, a full disk) must not wedge a
    gene forever: a decode/OS error logs a warning and returns None so the gene
    simply recomputes (and overwrites the bad file) instead of raising a
    non-retryable JSONDecodeError on every future run. Mirrors the guarded reads
    in src/pathways.py, src/fetchers/cellxgene.py, and scripts/build_synonyms.py.
    """
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("Cache file %s is unreadable (%s); ignoring it", path, exc)
        return None


def read_cache(gene: str, config: PipelineConfig) -> dict | None:
    """Return a gene's final enriched record, or None on miss / disabled / force.

    Honors ENABLE_CACHE (off → no cache), FORCE_RERUN (recompute everything), and
    FORCE_REMERGE (recompute merge + enrich from the raw cache). All three cause a
    final-cache miss; the fresh result is still written back afterward.
    """
    if not config.enable_cache or config.force_rerun or config.force_remerge:
        return None
    path = _final_cache_path(gene, config)
    if path.exists():
        data = _load_cache_file(path)
        if data is not None:
            log.info("Gene %s: final-cache hit → skipping all stages", gene)
            return data
    return None


def write_cache(gene: str, config: PipelineConfig, result: dict) -> None:
    """Persist a gene's final enriched record (no-op if caching disabled)."""
    if not config.enable_cache:
        return
    path = _final_cache_path(gene, config)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, default=str))
    log.info("Gene %s: final result cached to %s", gene, path)


def read_raw_cache(gene: str, config: PipelineConfig) -> list[dict] | None:
    """Return a gene's cached per-source extractions (unfiltered), or None on miss.

    Honors ENABLE_CACHE and FORCE_RERUN (which bypasses both layers). FORCE_REMERGE
    deliberately does NOT bypass this layer — that is what lets a re-merge skip
    re-fetching/re-extracting: the extractions are reused and only the confidence
    filter + merge + enrich rerun. An empty list is a valid hit (the gene yielded
    no extraction at all) and is returned as ``[]``.
    """
    if not config.enable_cache or config.force_rerun:
        return None
    path = _raw_cache_path(gene, config)
    if path.exists():
        data = _load_cache_file(path)
        if data is not None:
            log.info("Gene %s: raw extraction cache hit", gene)
            return data.get("extractions")
    return None


def write_raw_cache(
    gene: str, config: PipelineConfig, extractions: list[dict]
) -> None:
    """Persist a gene's unfiltered per-source extractions (no-op if disabled)."""
    if not config.enable_cache:
        return
    path = _raw_cache_path(gene, config)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"extractions": extractions}, indent=2, default=str))
    log.info("Gene %s: raw extraction cached to %s", gene, path)


async def run_gene(
    gene: str, config: PipelineConfig, stats: RunStats
) -> dict | None:
    """Run all stages for a single gene through the two-layer resume cache.

    Execution order:
      1. Final-cache (full_key) hit → return the cached record, no work done
         (counts in ``stats.genes_cached``).
      2. Raw-cache (extract_key) hit → skip fetch + extract and replay merge +
         enrich from the cached extractions, no fetch/extract API cost (counts in
         ``stats.genes_remerged``). This is the fetch/extract-free path a synonym
         edit or FORCE_REMERGE takes (the merge model still runs for multi-source
         genes).
      3. Miss on both → run fetch + extract, write the raw cache, then merge +
         enrich.

    Returns the enriched annotation, or ``None`` when no high-confidence source
    survived the merge gate.
    """
    # Layer 2 — final enriched cache: a hit short-circuits the whole chain.
    cached = read_cache(gene, config)
    if cached is not None:
        stats.genes_cached += 1
        # Marker so the progress bar can show "cached"; stripped in main() before
        # the record is written to final_annotations.json.
        cached["_from_cache"] = True
        return cached

    # Layer 1 — raw extraction cache: a hit replays merge + enrich with no
    # fetch/extract API calls (the fetch/extract-free re-merge path).
    extractions = read_raw_cache(gene, config)
    if extractions is not None:
        log.info(
            "Gene %s: raw-cache hit → replaying merge + enrich "
            "(no fetch/extract API calls)",
            gene,
        )
        stats.genes_remerged += 1
    else:
        fetched = await run_fetch_stage(gene, config)
        extractions, extract_usages = await run_extract_stage(gene, fetched, config)
        for usage in extract_usages:
            stats.add_usage(usage)
        write_raw_cache(gene, config, extractions)

    # Apply the confidence gate here (not at extract time) so the raw cache holds
    # the unfiltered extractions: a CONFIDENCE_THRESHOLD change invalidates the
    # final layer (via full_key) and is honored by re-filtering on a re-merge,
    # without re-extracting.
    sources = [
        a for a in extractions
        if a.get("confidence", 0) >= config.confidence_threshold
    ]

    merged, merge_usage = await run_merge_stage(gene, sources, config)
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

    # Pre-two-layer caches were flat files directly under cache_dir; they are not
    # read by the new raw/ + final/ layers. Warn so they can be cleaned up; this
    # run repopulates both layers as genes recompute.
    legacy = [p for p in Path(config.cache_dir).glob("*.json") if p.is_file()]
    if legacy:
        log.warning(
            "Found %d legacy flat cache file(s) in %s from before the two-layer "
            "cache; they are ignored. Delete them once raw/ and final/ are "
            "populated by this run.",
            len(legacy), config.cache_dir,
        )

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
