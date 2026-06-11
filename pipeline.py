"""Main orchestrator (CLAUDE.md Step 8).

Fetches PubMed / UniProt / OpenTargets for each gene, extracts structured
annotations, merges them, then builds the target network and prioritized table.

Run:  python pipeline.py
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
from pathlib import Path

import httpx
from dotenv import load_dotenv

load_dotenv()

from src.utils import setup_logging, load_gene_list, load_ref_set
from src.fetchers.pubmed import search_pmids, fetch_abstracts
from src.fetchers.uniprot import fetch_uniprot
from src.fetchers.opentargets import fetch_opentargets
from src.extractor import (
    extract_from_text,
    extract_from_uniprot,
    extract_from_opentargets,
)
from src.merger import merge_annotations
from src.network import (
    build_target_network,
    compute_priority_scores,
    save_network,
    save_prioritized_tsv,
)

CONFIDENCE_THRESHOLD = float(os.getenv("CONFIDENCE_THRESHOLD", "0.65"))
DISEASE_FILTER = os.getenv("DISEASE_FILTER", "cancer")

ANNOTATIONS_JSONL = Path("outputs/annotations.jsonl")

log = logging.getLogger("bio_annot.pipeline")


async def process_gene(gene: str, reactome_ref: set, session) -> dict | None:
    """Fetch all sources, extract annotations, merge, return merged dict or None."""
    # 1. PubMed abstracts → extraction
    pmids = await search_pmids(gene)
    abstracts = await fetch_abstracts(pmids)
    pubmed_text = "\n\n".join(
        f"PMID:{a['pmid']}\n{a['abstract']}" for a in abstracts
    )
    pubmed_ann = extract_from_text(gene, pubmed_text, pmids)

    # 2. UniProt → extraction
    up_data = await fetch_uniprot(gene)
    uniprot_ann = extract_from_uniprot(gene, up_data) if up_data else None

    # 3. OpenTargets → extraction
    ot_data = await fetch_opentargets(gene)
    ot_ann = extract_from_opentargets(gene, ot_data) if ot_data else None

    # 4. Keep only high-confidence sources
    sources = [
        a
        for a in [pubmed_ann, uniprot_ann, ot_ann]
        if a and a.get("confidence", 0) >= CONFIDENCE_THRESHOLD
    ]

    # 5. Bail if nothing survives the confidence gate
    if not sources:
        log.warning("No high-confidence sources for %s", gene)
        return None

    # 6. Merge
    merged = merge_annotations(gene, sources, reactome_ref)

    # 7. Persist raw sources and append the merged record
    raw_path = Path("outputs/raw") / f"{gene}_raw.json"
    with open(raw_path, "w", encoding="utf-8") as f:
        json.dump(
            {"pubmed": pubmed_ann, "uniprot": uniprot_ann, "opentargets": ot_ann},
            f,
            indent=2,
        )
    with open(ANNOTATIONS_JSONL, "a", encoding="utf-8") as f:
        f.write(json.dumps(merged) + "\n")

    # 8. Done
    return merged


async def main() -> None:
    setup_logging(os.getenv("LOG_LEVEL", "INFO"))
    Path("outputs/raw").mkdir(parents=True, exist_ok=True)
    # Start each run with a fresh annotations log so re-runs don't accumulate.
    ANNOTATIONS_JSONL.unlink(missing_ok=True)

    genes = load_gene_list("inputs/target_genes.txt")
    reactome_ref = load_ref_set("refs/reactome_pathways.txt")
    log.info("Processing %d genes | disease_filter=%s", len(genes), DISEASE_FILTER)

    final_annotations: dict = {}
    # Process genes with a concurrency limit of 3 (respect API rate limits).
    semaphore = asyncio.Semaphore(3)

    async def bounded(gene: str) -> None:
        async with semaphore:
            result = await process_gene(gene, reactome_ref, session=None)
            if result:
                final_annotations[gene] = result

    async with httpx.AsyncClient(timeout=30.0) as session:
        await asyncio.gather(*[bounded(g) for g in genes])

    # Write final merged JSON
    with open("outputs/final_annotations.json", "w", encoding="utf-8") as f:
        json.dump(final_annotations, f, indent=2)
    log.info(
        "Wrote %d annotations → outputs/final_annotations.json", len(final_annotations)
    )

    # Build network and prioritize
    G = build_target_network(final_annotations)
    save_network(G, "outputs/target_network.gpickle")
    scores = compute_priority_scores(G, DISEASE_FILTER)
    save_prioritized_tsv(scores, "outputs/prioritized_targets.tsv")
    log.info("Top 5 targets: %s", [s["gene"] for s in scores[:5]])


if __name__ == "__main__":
    asyncio.run(main())
