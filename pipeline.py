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
from pathlib import Path
from typing import Any

from dotenv import load_dotenv

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
) -> list[dict]:
    """Extract structured annotations from each source.

    Runs the LLM extractor over PubMed text, UniProt, and OpenTargets, persists
    the (unfiltered) per-source extractions to ``outputs/raw/{gene}_raw.json``,
    and returns only the sources clearing the confidence gate.
    """
    pubmed_ann = await extract_from_text(
        gene, fetched["pubmed_text"], fetched["pmids"]
    )
    up_data = fetched["uniprot"]
    ot_data = fetched["opentargets"]
    uniprot_ann = await extract_from_uniprot(gene, up_data) if up_data else None
    ot_ann = await extract_from_opentargets(gene, ot_data) if ot_data else None

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
    return sources


async def run_merge_stage(
    gene: str, extractions: list[dict], config: PipelineConfig
) -> dict | None:
    """Merge and resolve conflicts across sources.

    Returns ``None`` (a clean skip, logged) when no source cleared the
    confidence gate, so the gene is excluded from output without being treated
    as a failure.
    """
    if not extractions:
        log.warning("No high-confidence sources for %s", gene)
        return None
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


async def run_gene(gene: str, config: PipelineConfig) -> dict | None:
    """Run all stages for a single gene.

    Returns the enriched annotation, or ``None`` when no high-confidence source
    survived the merge gate. On a resume-cache hit, returns the cached result
    without running any stage.
    """
    # Check the on-disk resume cache first.
    cached = read_cache(gene, config)
    if cached is not None:
        return cached

    fetched = await run_fetch_stage(gene, config)
    extractions = await run_extract_stage(gene, fetched, config)
    merged = await run_merge_stage(gene, extractions, config)
    if merged is None:
        return None
    enriched = await run_enrich_stage(gene, merged, config)

    # Append the merged record to the run's annotations log, then cache it.
    with open(ANNOTATIONS_JSONL, "a", encoding="utf-8") as f:
        f.write(json.dumps(enriched) + "\n")
    write_cache(gene, config, enriched)
    return enriched


async def main() -> None:
    setup_logging(config.log_level)
    Path("outputs/raw").mkdir(parents=True, exist_ok=True)
    # Start each run with a fresh annotations log so re-runs don't accumulate.
    ANNOTATIONS_JSONL.unlink(missing_ok=True)

    log.info("Running in disease context: %s", config.disease_context)
    genes = load_gene_list("inputs/target_genes.txt")
    log.info("Processing %d genes", len(genes))

    # Process genes with a concurrency limit (respect API rate limits).
    sem = asyncio.Semaphore(config.semaphore_limit)

    async def bounded(gene: str) -> dict | None:
        async with sem:
            return await run_gene(gene, config)

    results = await asyncio.gather(
        *[bounded(g) for g in genes], return_exceptions=True
    )

    # Error isolation: one gene's failure (e.g. a ConnectTimeout) must not kill
    # the whole run — log it, exclude it, and keep the rest.
    final_annotations: dict = {}
    failed = []
    for gene, result in zip(genes, results):
        if isinstance(result, Exception):
            log.error("Gene %s failed: %s, skipping", gene, result)
            failed.append(gene)
        elif result:
            final_annotations[gene] = result
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


if __name__ == "__main__":
    asyncio.run(main())
