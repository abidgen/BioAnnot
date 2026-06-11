"""Batch API variant for >= 50 genes (CLAUDE.md Step 9).

Fetches all sources synchronously, submits one Anthropic Message Batch
(~50% cheaper), polls until it ends, then merges (single-source) and builds
the same network/prioritization outputs as pipeline.py.

Run:  python batch_pipeline.py
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import time
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

from src.utils import setup_logging, load_gene_list, load_ref_set
from src.fetchers.pubmed import search_pmids, fetch_abstracts
from src.fetchers.uniprot import fetch_uniprot
from src.fetchers.opentargets import fetch_opentargets
from src.extractor import (
    ANNOTATION_TOOL,
    SYSTEM_PROMPT,
    _get_client,
    _truncate_words,
    _format_uniprot_text,
    _format_opentargets_text,
)
from src.merger import merge_annotations
from src.network import (
    build_target_network,
    compute_priority_scores,
    save_network,
    save_prioritized_tsv,
)

BATCH_MODEL = "claude-opus-4-8"
MAX_TOKENS = 1024
POLL_INTERVAL_SECONDS = 60
DISEASE_FILTER = os.getenv("DISEASE_FILTER", "cancer")

log = logging.getLogger("bio_annot.batch")


async def _fetch_gene(gene: str, semaphore: asyncio.Semaphore) -> tuple[str, str]:
    """Fetch all sources for one gene and assemble its combined input text."""
    async with semaphore:
        pmids = await search_pmids(gene)
        abstracts = await fetch_abstracts(pmids)
        pubmed_block = "\n\n".join(
            f"PMID:{a['pmid']}\n{a['abstract']}" for a in abstracts
        )

        up_data = await fetch_uniprot(gene)
        uniprot_text = _format_uniprot_text(up_data) if up_data else ""

        ot_data = await fetch_opentargets(gene)
        ot_text = _format_opentargets_text(ot_data) if ot_data else ""

    combined_text = "\n\n".join([pubmed_block, uniprot_text, ot_text])
    # Truncate to the same ~3000-word budget the extractor enforces, so a
    # gene with lots of literature can't overflow the batch request context.
    return gene, _truncate_words(combined_text)


async def fetch_all(genes: list[str]) -> list[tuple[str, str]]:
    """Fetch (gene, combined_text) pairs for every gene, 3 at a time."""
    semaphore = asyncio.Semaphore(3)
    return await asyncio.gather(*[_fetch_gene(g, semaphore) for g in genes])


def main() -> None:
    setup_logging(os.getenv("LOG_LEVEL", "INFO"))
    Path("outputs/raw").mkdir(parents=True, exist_ok=True)

    genes = load_gene_list("inputs/target_genes.txt")
    reactome_ref = load_ref_set("refs/reactome_pathways.txt")
    log.info("Batch processing %d genes | disease_filter=%s", len(genes), DISEASE_FILTER)

    # 1. Fetch all source text synchronously.
    gene_texts = asyncio.run(fetch_all(genes))

    client = _get_client()

    # 2. Build batch requests.
    requests = [
        {
            "custom_id": f"gene-{gene}",
            "params": {
                "model": BATCH_MODEL,
                "max_tokens": MAX_TOKENS,
                "tools": ANNOTATION_TOOL,
                "tool_choice": {"type": "tool", "name": "annotate_target"},
                "system": SYSTEM_PROMPT,
                "messages": [
                    {"role": "user", "content": f"Gene: {gene}\n\n{combined_text}"}
                ],
            },
        }
        for gene, combined_text in gene_texts
    ]

    # 3. Submit and persist the batch id.
    batch = client.messages.batches.create(requests=requests)
    log.info("Submitted batch %s", batch.id)
    print(f"Batch ID: {batch.id}")
    with open("outputs/batch_id.txt", "w", encoding="utf-8") as f:
        f.write(batch.id + "\n")

    # 4. Poll until the batch ends.
    while batch.processing_status != "ended":
        time.sleep(POLL_INTERVAL_SECONDS)
        batch = client.messages.batches.retrieve(batch.id)
        log.info("Batch %s status=%s counts=%s", batch.id, batch.processing_status, batch.request_counts)

    # 5. Collect results.
    raw_annotations: dict = {}
    for result in client.messages.batches.results(batch.id):
        gene = result.custom_id.replace("gene-", "")
        if result.result.type != "succeeded":
            log.warning("Batch request for %s did not succeed: %s", gene, result.result.type)
            continue
        tool_use = next(
            (b for b in result.result.message.content if b.type == "tool_use"), None
        )
        if tool_use is None:
            log.warning("No tool_use block in result for %s", gene)
            continue
        raw_annotations[gene] = [tool_use.input]  # single source; no merge needed

    # 6. Merge pass (single-source merge just validates pathway names).
    final_annotations: dict = {}
    for gene, sources in raw_annotations.items():
        merged = merge_annotations(gene, sources, reactome_ref)
        final_annotations[gene] = merged
        with open(Path("outputs/raw") / f"{gene}_raw.json", "w", encoding="utf-8") as f:
            json.dump(sources, f, indent=2)

    with open("outputs/final_annotations.json", "w", encoding="utf-8") as f:
        json.dump(final_annotations, f, indent=2)
    log.info(
        "Wrote %d annotations → outputs/final_annotations.json", len(final_annotations)
    )

    # 7. Build network and prioritize — same as pipeline.py.
    G = build_target_network(final_annotations)
    save_network(G, "outputs/target_network.gpickle")
    scores = compute_priority_scores(G, DISEASE_FILTER)
    save_prioritized_tsv(scores, "outputs/prioritized_targets.tsv")
    log.info("Top 5 targets: %s", [s["gene"] for s in scores[:5]])


if __name__ == "__main__":
    main()
