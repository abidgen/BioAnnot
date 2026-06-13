"""Multi-source merge with LLM-assisted conflict resolution (CLAUDE.md Step 6).

Uses claude-sonnet-4-6 (cheaper than the extraction model) with the same
annotate_target tool to reconcile per-source annotations into one record, then
validates pathway names against the canonical Reactome reference set.
"""

from __future__ import annotations

import json
import logging
import os
import re
from datetime import datetime, timezone
from pathlib import Path

from rapidfuzz import fuzz, process

from src.extractor import ANNOTATION_TOOL, _get_client
from src.utils import CACHE_MIN_TOKENS, estimate_tokens, load_ref_set, validate_pmids

log = logging.getLogger("bio_annot.merger")

# Model is env-configurable; default preserves the CLAUDE.md merge model.
MERGE_MODEL = os.getenv("MERGE_MODEL", "claude-sonnet-4-6")
MAX_TOKENS = 4096

NON_CANONICAL_PREFIX = "NON-CANONICAL: "

# Minimum rapidfuzz token_sort_ratio (0–100) for a pathway name to be accepted
# as a fuzzy match to a canonical Reactome name.
FUZZY_THRESHOLD = float(os.getenv("FUZZY_THRESHOLD", "85"))

# Informal→canonical pathway synonym map (built by scripts/build_synonyms.py) and
# the canonical Reactome reference used to validate it. The map is loaded and
# validated below, once _normalize is defined.
SYNONYMS_PATH = Path("refs/pathway_synonyms.json")
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

# Reactome stable-ID suffix the merge model often appends to a pathway name,
# e.g. "Oxidative Stress Induced Senescence (R-HSA-2559580)". The OpenTargets
# extractor (extractor.py _format_opentargets_text) feeds names in this form, so
# the model copies it. The canonical reference stores bare names, so we strip
# the suffix before comparing.
_RHSA_SUFFIX = re.compile(r"\s*\(R-HSA-\d+\)\s*$")


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


def _normalize(name: str) -> str:
    """Normalize a pathway name for canonical comparison.

    Lowercases, trims, and strips any trailing Reactome stable-ID suffix so a
    name decorated with "(R-HSA-…)" still matches its bare canonical form.
    """
    return _RHSA_SUFFIX.sub("", name.strip()).strip().lower()


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


def _load_synonyms(reactome_ref: set[str]) -> dict[str, str]:
    """Load and validate the informal→canonical pathway synonym map.

    Keys are lowercased and null/empty values skipped. Each non-null value is
    validated against ``reactome_ref`` (normalized comparison): a value that is
    not a real Reactome name is dropped with a WARNING and treated as null — so a
    synonym lookup can only ever yield a genuine canonical name, returned in the
    reference's exact casing. This guards against the synonym-builder LLM
    emitting plausible-sounding but non-existent Reactome names. If the Reactome
    reference is unavailable, validation is skipped (values loaded as-is).
    """
    if not SYNONYMS_PATH.exists():
        return {}
    try:
        with open(SYNONYMS_PATH, encoding="utf-8") as f:
            raw = json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("Could not load %s (%s); fuzzy synonyms disabled", SYNONYMS_PATH, exc)
        return {}

    if not reactome_ref:
        log.warning(
            "Reactome reference unavailable; loading %d synonym(s) without validation",
            sum(1 for v in raw.values() if v),
        )
        return {k.lower(): v for k, v in raw.items() if v}

    canonical_by_norm = {_normalize(r): r for r in reactome_ref}
    validated: dict[str, str] = {}
    for informal, value in raw.items():
        if not value:
            continue
        canonical = canonical_by_norm.get(_normalize(value))
        if canonical is None:
            log.warning(
                "Synonym %r → %r is not in the Reactome reference; ignoring it",
                informal,
                value,
            )
            continue
        validated[informal.lower()] = canonical
    return validated


_REACTOME_REF = _load_reactome_ref()
_SYNONYMS = _load_synonyms(_REACTOME_REF)


def _utc_now_iso() -> str:
    """Current UTC time as an ISO-8601 string with a trailing Z."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# Gene/disease-templated Reactome pathway names. The gene (or gene-pair) symbol
# is only a few characters of a long shared template, so a fuzzy scorer rates
# wrong-gene siblings (e.g. "Signaling by BRCA1 mutants" vs "Signaling by AMER1
# mutants") highly. Each entry is (regex, capture-group-index); the captured
# token is the gene/key symbol used by the guard. Order matters — the first
# matching pattern wins. Kept in sync with the copy in scripts/build_synonyms.py
# (that standalone script can't import this module).
TEMPLATED_PATTERNS = [
    # Signaling by <GENE> (including ligand-responsive variants)
    (r'(?i)^(?:constitutive\s+)?signaling by (?:ligand-responsive\s+)?([A-Z][\w\d]*/[\w\d]+)', 1),
    (r'(?i)^(?:constitutive\s+)?signaling by (?:ligand-responsive\s+)?([A-Z][\w\d]*)', 1),

    # Signaling by <GENE> in cancer / mutants
    (r'(?i)^signaling by ([A-Z][\w\d]*)\s+(in cancer|mutants?)', 1),

    # Nuclear events stimulated by <GENE>
    (r'(?i)^nuclear events stimulated by ([A-Z][\w\d]*)', 1),

    # Defective <GENE> causes <DISEASE>
    (r'(?i)^defective ([A-Z][\w\d]*)\s+causes', 1),

    # <GENE> variants cause <DISEASE>
    (r'(?i)^([A-Z][\w\d]*)\s+variants?\s+cause', 1),

    # Loss of Function of <GENE>
    (r'(?i)^loss of (?:function of\s+)?([A-Z][\w\d]*/[\w\d]+)', 1),
    (r'(?i)^loss of (?:function of\s+)?([A-Z][\w\d]*)', 1),

    # <GENE> Loss of Function in Cancer
    (r'(?i)^([A-Z][\w\d]*)\s+loss of function in cancer', 1),

    # Regulation of <GENE> activity/signaling/expression/function/degradation
    (r'(?i)^regulation of ([A-Z][\w\d]*)\s+(activity|signaling|expression|degradation|function)', 1),

    # Activation of <GENE>
    (r'(?i)^activation of ([A-Z][\w\d]*)', 1),

    # <GENE> mediated ...
    (r'(?i)^([A-Z][\w\d]*)\s+mediated', 1),

    # Drug resistance of/in <GENE> mutants
    (r'(?i)^drug resistance (?:of|in) ([A-Z][\w\d]*)', 1),

    # <DRUG>-resistant <GENE> mutants
    (r'(?i)^\w+-resistant ([A-Z][\w\d]*)\s+mutants?', 1),

    # Slash-separated gene pairs: SMAD2/3, PI3K/AKT
    (r'(?i)^([A-Z][\w\d]*/[A-Z][\w\d]*)\s+', 1),

    # Aberrant regulation of ... due to <GENE> defects
    (r'(?i)^aberrant regulation of .+ due to ([A-Z][\w\d]*)\s+defects', 1),
]


def extract_key_token(name: str) -> str | None:
    """Return the gene/key token (uppercased) from a templated name, or None."""
    for pattern, group in TEMPLATED_PATTERNS:
        m = re.match(pattern, name.strip())
        if m:
            return m.group(group).upper()
    return None


def gene_token_guard(query: str, candidate: str) -> bool:
    """Reject a fuzzy match whose key token differs from the query's.

    Returns False only when both names are templated and their tokens differ;
    otherwise (one/neither templated, or tokens equal) the match is allowed.
    """
    q_token = extract_key_token(query)
    c_token = extract_key_token(candidate)
    if q_token and c_token and q_token != c_token:
        return False
    return True


def _fuzzy_canonical(
    pathway: str, reactome_ref: set[str]
) -> tuple[bool, str, str]:
    """Resolve a pathway name to its canonical Reactome form.

    Returns ``(is_canonical, method, name)`` where ``method`` is one of
    ``exact`` / ``synonym`` / ``fuzzy`` / ``non_canonical``. On a canonical
    match ``name`` is the exact Reactome string; otherwise it is the bare input.
    Resolution proceeds in priority order:

      a) strip the R-HSA suffix, then a normalized exact match → "exact";
      b) the informal→canonical synonym map (case-insensitive) → "synonym";
      c) rapidfuzz token_sort_ratio ≥ FUZZY_THRESHOLD → "fuzzy";
      d) otherwise → "non_canonical".

    Logs the method (and score) for every call.
    """
    bare = (
        pathway[len(NON_CANONICAL_PREFIX):]
        if pathway.startswith(NON_CANONICAL_PREFIX)
        else pathway
    ).strip()
    key = _normalize(bare)

    # a) Exact match (normalized) → return the reference's canonical casing.
    canonical_by_norm = {_normalize(r): r for r in reactome_ref}
    if key in canonical_by_norm:
        canonical = canonical_by_norm[key]
        log.info("pathway %r → exact %r (score=100.0)", bare, canonical)
        return True, "exact", canonical

    # b) Synonym map (lowercased keys; null entries already dropped on load).
    synonym = _SYNONYMS.get(bare.lower())
    if synonym:
        log.info("pathway %r → synonym %r", bare, synonym)
        return True, "synonym", synonym

    # c) Fuzzy match against the full canonical set (case-insensitive scoring).
    choices = list(reactome_ref)
    best = (
        process.extractOne(
            bare, choices, scorer=fuzz.token_sort_ratio, processor=str.lower
        )
        if choices
        else None
    )
    if best is not None:
        match, score, _ = best
        if score >= FUZZY_THRESHOLD and gene_token_guard(bare, match):
            log.info("pathway %r → fuzzy %r (score=%.1f)", bare, match, score)
            return True, "fuzzy", match
        if score >= FUZZY_THRESHOLD:
            # Score cleared but the gene token differs — reject the sibling. With
            # no LLM at merge time, this falls through to non_canonical.
            log.info(
                "pathway %r → non_canonical (fuzzy %r score=%.1f rejected by "
                "gene-token guard)",
                bare,
                match,
                score,
            )
        else:
            log.info(
                "pathway %r → non_canonical (best fuzzy %r score=%.1f < %.0f)",
                bare,
                match,
                score,
                FUZZY_THRESHOLD,
            )
    else:
        log.info("pathway %r → non_canonical (no reactome reference)", bare)

    # d) No canonical resolution.
    return False, "non_canonical", bare


def _flag_noncanonical(pathways: list[str], reactome_ref: set[str]) -> list[str]:
    """Canonicalize pathway names; prefix unresolved ones with NON-CANONICAL:.

    Each name is run through :func:`_fuzzy_canonical`: exact/synonym/fuzzy hits
    are replaced with the exact Reactome string, everything else is prefixed
    ``NON-CANONICAL: ``. Idempotent — a name already carrying the prefix is
    checked on its bare form, never double-prefixed.
    """
    flagged: list[str] = []
    for pathway in pathways:
        is_canonical, _method, name = _fuzzy_canonical(pathway, reactome_ref)
        flagged.append(name if is_canonical else NON_CANONICAL_PREFIX + name)
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

    response = _get_client().messages.create(
        model=MERGE_MODEL,
        max_tokens=MAX_TOKENS,
        system=[
            {
                "type": "text",
                "text": MERGE_SYSTEM_PROMPT,
                "cache_control": {"type": "ephemeral"},
            }
        ],
        tools=ANNOTATION_TOOL,
        tool_choice={"type": "tool", "name": "annotate_target"},
        messages=[
            {
                "role": "user",
                "content": f"Gene: {gene} ({n} sources)\n\n{sources_block}",
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
    log.info(
        "merge_annotations(%s) cache: %s read, %s created",
        gene,
        usage.cache_read_input_tokens,
        usage.cache_creation_input_tokens,
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
