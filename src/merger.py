"""Multi-source merge with LLM-assisted conflict resolution (CLAUDE.md Step 6).

Uses claude-sonnet-4-6 (cheaper than the extraction model) with the same
annotate_target tool to reconcile per-source annotations into one record, then
validates pathway names against the canonical Reactome reference set.
"""

from __future__ import annotations

import json
import logging
import random
import re
from datetime import datetime, timezone

from src.extractor import ANNOTATION_TOOL, _get_client
from src.utils import validate_pmids

log = logging.getLogger("bio_annot.merger")

MERGE_MODEL = "claude-sonnet-4-6"
MAX_TOKENS = 4096

NON_CANONICAL_PREFIX = "NON-CANONICAL: "

# Number of valid Reactome names to show the model as exact-naming examples.
PATHWAY_EXAMPLE_COUNT = 20

# Reactome stable-ID suffix the merge model often appends to a pathway name,
# e.g. "Oxidative Stress Induced Senescence (R-HSA-2559580)". The OpenTargets
# extractor (extractor.py _format_opentargets_text) feeds names in this form, so
# the model copies it. The canonical reference stores bare names, so we strip
# the suffix before comparing.
_RHSA_SUFFIX = re.compile(r"\s*\(R-HSA-\d+\)\s*$")


def _normalize(name: str) -> str:
    """Normalize a pathway name for canonical comparison.

    Lowercases, trims, and strips any trailing Reactome stable-ID suffix so a
    name decorated with "(R-HSA-…)" still matches its bare canonical form.
    """
    return _RHSA_SUFFIX.sub("", name.strip()).strip().lower()


def _utc_now_iso() -> str:
    """Current UTC time as an ISO-8601 string with a trailing Z."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _flag_noncanonical(pathways: list[str], reactome_ref: set[str]) -> list[str]:
    """Prefix any pathway not in the canonical Reactome set with NON-CANONICAL:.

    The match is case- and whitespace-insensitive (both sides normalized), so
    "EGFR Signaling Pathway" matches "EGFR signaling pathway". The original
    casing of the model's name is preserved in the output. Idempotent — a name
    already carrying the prefix is checked on its bare form, never double-prefixed.
    """
    normalized_ref = {_normalize(r) for r in reactome_ref}
    flagged: list[str] = []
    for pathway in pathways:
        bare = (
            pathway[len(NON_CANONICAL_PREFIX):]
            if pathway.startswith(NON_CANONICAL_PREFIX)
            else pathway
        )
        if _normalize(bare) in normalized_ref:
            flagged.append(bare)
        else:
            flagged.append(NON_CANONICAL_PREFIX + bare)
    return flagged


def _union_source_pmids(source_annotations: list[dict]) -> list[str]:
    """Union (validated, deduped, sorted) source_pmids across all sources."""
    pmids: set[str] = set()
    for source in source_annotations:
        pmids.update(validate_pmids(source.get("source_pmids", [])))
    return sorted(pmids)


def merge_annotations(
    gene: str, source_annotations: list[dict], reactome_ref: set[str]
) -> dict:
    """Reconcile per-source annotation records for a gene into one merged record.

    With a single source, returns it directly (validating pathway names). With
    multiple sources, calls the merge model with the annotate_target tool forced.
    """
    n = len(source_annotations)

    if n == 1:
        merged = dict(source_annotations[0])
        merged["pathways"] = _flag_noncanonical(
            merged.get("pathways", []), reactome_ref
        )
        merged["source_pmids"] = _union_source_pmids(source_annotations)
        merged["source_count"] = 1
        merged["merged_at"] = _utc_now_iso()
        return merged

    sources_block = json.dumps(source_annotations, indent=2)

    # Show the model a deterministic sample of real Reactome names so it emits
    # exact canonical names (seeded for prompt-cache stability).
    ref_sorted = sorted(reactome_ref)
    sample = random.Random(0).sample(
        ref_sorted, min(PATHWAY_EXAMPLE_COUNT, len(ref_sorted))
    )
    pathway_examples = "; ".join(sample)

    system_prompt = (
        f"You are reconciling annotation records for {gene} from {n} sources.\n"
        "When naming pathways, you MUST use exact Reactome pathway names. "
        f"Valid examples include: {pathway_examples}\n"
        "RULES:\n"
        "- Include a pathway only if ≥2 sources agree OR a single source rates "
        "confidence ≥ 0.85\n"
        "- For conflicting disease roles (e.g. oncogene vs suppressor), note "
        "context-dependence in the role field\n"
        "- Union all interactors from all sources\n"
        "- Set final confidence = mean(source confidences) × "
        "(1 - 0.1 × conflict_count)\n"
        "- Flag any pathway name not found in the canonical Reactome list with "
        f"prefix '{NON_CANONICAL_PREFIX}'\n"
        "- Only assert what is supported by at least one source"
    )

    response = _get_client().messages.create(
        model=MERGE_MODEL,
        max_tokens=MAX_TOKENS,
        system=system_prompt,
        tools=ANNOTATION_TOOL,
        tool_choice={"type": "tool", "name": "annotate_target"},
        messages=[
            {
                "role": "user",
                "content": f"Source annotations for {gene}:\n\n{sources_block}",
            }
        ],
    )

    usage = response.usage
    log.info(
        "merge_annotations(%s): %d sources, input_tokens=%s output_tokens=%s",
        gene,
        n,
        usage.input_tokens,
        usage.output_tokens,
    )

    tool_use = next(
        (block for block in response.content if block.type == "tool_use"), None
    )
    if tool_use is None:
        log.warning("No tool_use block returned when merging %s", gene)
        merged = {"gene_symbol": gene, "functions": [], "pathways": [], "confidence": 0.0}
    else:
        merged = dict(tool_use.input)

    # Validate pathway names against the canonical Reactome reference.
    merged["pathways"] = _flag_noncanonical(merged.get("pathways", []), reactome_ref)
    # Propagate PMIDs: the merge tool schema has no source_pmids field, so union
    # them from the per-source records rather than relying on the model output.
    merged["source_pmids"] = _union_source_pmids(source_annotations)
    merged["source_count"] = n
    merged["merged_at"] = _utc_now_iso()
    return merged
