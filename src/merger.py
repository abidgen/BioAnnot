"""Multi-source merge with LLM-assisted conflict resolution (CLAUDE.md Step 6).

Uses claude-sonnet-4-6 (cheaper than the extraction model) with the same
annotate_target tool to reconcile per-source annotations into one record, then
validates pathway names against the canonical Reactome reference set.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from typing import Any

from src.config import config
from src.extractor import ANNOTATION_TOOL
from src.llm import call_tool
from src.pathways import NON_CANONICAL_PREFIX, fuzzy_canonical, load_synonyms
from src.utils import CACHE_MIN_TOKENS, estimate_tokens, load_ref_set, validate_pmids

log = logging.getLogger("bio_annot.merger")

# Model and token budget are centralized in src.config (env-configurable).
MERGE_MODEL = config.merge_model
MAX_TOKENS = config.max_tokens

# Minimum rapidfuzz token_sort_ratio (0–100) for a pathway name to be accepted
# as a fuzzy match to a canonical Reactome name.
FUZZY_THRESHOLD = config.fuzzy_threshold

# Canonical Reactome reference used to validate the synonym map at import.
REACTOME_PATH = "refs/reactome_pathways.txt"

# Static merge system prompt. Per-gene context (gene symbol, source count, source
# records) lives entirely in the user message, so this prompt is identical on
# every call and forms a cacheable prefix together with the tool schema.
MERGE_SYSTEM_PROMPT = (
    "You are a biological annotation merge specialist. You reconcile annotation "
    "records from multiple sources into a single high-quality structured "
    "annotation, resolving conflicts by evidence strength and source agreement.\n\n"
    "You will receive several per-source annotation records for one gene (each a "
    "JSON object derived from PubMed, UniProt, or OpenTargets); the user message "
    "names the gene and the number of sources. You must emit exactly one merged "
    "record via the annotate_target tool.\n\n"
    "GENERAL RULES:\n"
    "- Only assert what is supported by at least one source.\n"
    "- Include a pathway only if two or more sources agree OR a single source "
    "rates confidence >= 0.85.\n"
    "- Union all interactors across sources; deduplicate by gene symbol.\n"
    "- Flag any pathway not found in the canonical Reactome list with the prefix "
    "'NON-CANONICAL: '.\n\n"
    "PATHWAY NAMING:\n"
    "Always use exact Reactome pathway names. Valid examples: Signaling by WNT, "
    "RAF/MAP kinase cascade, Transcriptional Regulation by TP53. Prefer the "
    "canonical name over informal shorthand whenever one exists.\n\n"
    "CONFLICT RESOLUTION BY FIELD:\n"
    "- functions: When sources disagree, prefer the more specific, mechanistically "
    "precise description. Source priority for function: UniProt > OpenTargets > "
    "PubMed. UniProt curation is authoritative for molecular function and "
    "subcellular localization; defer to it.\n"
    "- pathways: Source priority: Reactome > OpenTargets > PubMed. Prefer canonical "
    "Reactome names from the higher-priority source and drop informal pathway "
    "shorthand when a canonical name is available.\n"
    "- disease_associations: For conflicting disease roles (e.g. oncogene vs "
    "tumor_suppressor), do NOT pick one — note the context-dependence in the role "
    "field. Keep the highest evidence_strength supported across sources.\n"
    "- cellular_states: Union across sources; deduplicate semantically equivalent "
    "terms.\n"
    "- interactors: Union across all sources, gene symbols only.\n"
    "- druggability_notes: Concatenate complementary detail; prefer the source with "
    "concrete drug or binding-pocket information.\n\n"
    "CONFIDENCE MERGING:\n"
    "Set the final confidence to a weighted average of the source confidences, "
    "weighted by source reliability (curated UniProt and OpenTargets records weigh "
    "more than PubMed free-text extraction), then reduce it by 0.1 for each "
    "material conflict you had to resolve. Clamp the result to the range 0.0 to "
    "1.0. A material conflict is a direct disagreement on a substantive claim — a "
    "contradictory disease role, an incompatible molecular function, or a pathway "
    "asserted by one source and explicitly excluded by another — not a mere "
    "difference in how many items each source happened to list.\n\n"
    "DEDUPLICATION AND PHRASING:\n"
    "Merge semantically equivalent entries even when their wording differs (for "
    "example 'sequence-specific DNA-binding transcription factor' and 'transcription "
    "factor activity' describe the same function — keep the more precise one). "
    "Preserve canonical gene symbols (HGNC) and exact Reactome pathway names. Keep "
    "each list item a short phrase rather than a sentence, and cap functions at "
    "eight items, retaining the best-supported ones.\n\n"
    "OUTPUT REQUIREMENTS:\n"
    "Never introduce a fact that appears in none of the sources, and never drop a "
    "well-supported fact simply to shorten the record. The merged record must be "
    "strictly grounded in the supplied sources. Emit exactly one annotate_target "
    "tool call representing the reconciled record."
)


def _log_cache_prefix_size() -> None:
    """Confirm the static system+tools prefix clears the prompt-cache minimum.

    Logged once at import; an undersized prefix is a WARNING so it surfaces even
    before setup_logging configures handlers.
    """
    prefix = MERGE_SYSTEM_PROMPT + json.dumps(ANNOTATION_TOOL)
    sys_tokens = estimate_tokens(MERGE_SYSTEM_PROMPT)
    prefix_tokens = estimate_tokens(prefix)
    level = logging.INFO if prefix_tokens >= CACHE_MIN_TOKENS else logging.WARNING
    log.log(
        level,
        "merger system prompt ~%d tokens; system+tools cache prefix ~%d tokens "
        "(cache min %d)",
        sys_tokens,
        prefix_tokens,
        CACHE_MIN_TOKENS,
    )


_log_cache_prefix_size()


def _load_reactome_ref() -> set[str]:
    """Load the canonical Reactome name set used to validate synonyms at import.

    Best-effort: returns an empty set (validation then skipped) if the reference
    file is missing or unreadable, so importing the module never hard-fails.
    """
    try:
        return load_ref_set(REACTOME_PATH)
    except OSError as exc:
        log.warning("Could not load Reactome reference %s (%s)", REACTOME_PATH, exc)
        return set()


# Canonical Reactome reference and validated synonym map, loaded once at import.
# Pathway normalization, the synonym loader, and fuzzy resolution all live in
# src.pathways (shared with scripts/build_synonyms.py).
_REACTOME_REF = _load_reactome_ref()
_SYNONYMS = load_synonyms(_REACTOME_REF)


def _utc_now_iso() -> str:
    """Current UTC time as an ISO-8601 string with a trailing Z."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _flag_noncanonical(pathways: list[str], reactome_ref: set[str]) -> list[str]:
    """Canonicalize pathway names; prefix unresolved ones with NON-CANONICAL:.

    Each name is run through :func:`src.pathways.fuzzy_canonical`: exact/synonym/
    fuzzy hits are replaced with the exact Reactome string, everything else is
    prefixed ``NON-CANONICAL: ``. Idempotent — a name already carrying the prefix
    is checked on its bare form, never double-prefixed.
    """
    flagged: list[str] = []
    for pathway in pathways:
        is_canonical, name, _method = fuzzy_canonical(
            pathway, reactome_ref, _SYNONYMS, FUZZY_THRESHOLD
        )
        flagged.append(name if is_canonical else NON_CANONICAL_PREFIX + name)
    return flagged


def _union_source_pmids(source_annotations: list[dict]) -> list[str]:
    """Union (validated, deduped, sorted) source_pmids across all sources."""
    pmids: set[str] = set()
    for source in source_annotations:
        pmids.update(validate_pmids(source.get("source_pmids", [])))
    return sorted(pmids)


async def merge_annotations(
    gene: str, source_annotations: list[dict[str, Any]], reactome_ref: set[str]
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Reconcile per-source annotation records for a gene into one merged record.

    With a single source, returns it directly (validating pathway names) with no
    LLM call. With multiple sources, calls the merge model with the
    annotate_target tool forced. Returns ``(merged, usage)`` where ``usage`` is
    the merge call's token-usage dict, or ``None`` for the single-source path
    (no LLM call was made).
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
        return merged, None

    sources_block = json.dumps(source_annotations, indent=2)

    merged_input, usage = await call_tool(
        model=MERGE_MODEL,
        system_prompt=MERGE_SYSTEM_PROMPT,
        messages=[
            {
                "role": "user",
                "content": f"Gene: {gene} ({n} sources)\n\n{sources_block}",
            }
        ],
        tools=ANNOTATION_TOOL,
        tool_name="annotate_target",
        max_tokens=MAX_TOKENS,
        label=f"merge_annotations({gene}, {n} sources)",
        # Pin to 0 for run-to-run reproducibility: a re-merge (e.g. after a
        # cache-invalidating config change) must reconcile the same sources the
        # same way, so the merged pathway set / network doesn't drift.
        temperature=0.0,
    )
    merged = dict(merged_input)

    # Validate pathway names against the canonical Reactome reference.
    merged["pathways"] = _flag_noncanonical(merged.get("pathways", []), reactome_ref)
    # Propagate PMIDs: the merge tool schema has no source_pmids field, so union
    # them from the per-source records rather than relying on the model output.
    merged["source_pmids"] = _union_source_pmids(source_annotations)
    merged["source_count"] = n
    merged["merged_at"] = _utc_now_iso()
    return merged, usage
