"""Shared pathway canonicalization logic.

Single home for the templated-pattern list, name normalization, the gene-token
guard, synonym-map loading, and fuzzy canonical resolution. Consumed by
``src.merger`` (at merge time) and ``scripts/build_synonyms.py`` (offline synonym
building) so there is exactly one implementation rather than two copies kept in
sync by hand.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

from rapidfuzz import fuzz, process

log = logging.getLogger("bio_annot.pathways")

NON_CANONICAL_PREFIX = "NON-CANONICAL: "

# Informal→canonical pathway synonym map (built by scripts/build_synonyms.py).
SYNONYMS_PATH = Path("refs/pathway_synonyms.json")

# Reactome stable-ID suffix the merge model often appends to a pathway name,
# e.g. "Oxidative Stress Induced Senescence (R-HSA-2559580)". The OpenTargets
# extractor feeds names in this form, so the model copies it. The canonical
# reference stores bare names, so we strip the suffix before comparing.
_RHSA_SUFFIX = re.compile(r"\s*\(R-HSA-\d+\)\s*$")


def _normalize(name: str) -> str:
    """Normalize a pathway name for canonical comparison.

    Lowercases, trims, and strips any trailing Reactome stable-ID suffix so a
    name decorated with "(R-HSA-…)" still matches its bare canonical form.
    """
    return _RHSA_SUFFIX.sub("", name.strip()).strip().lower()


# Gene/disease-templated Reactome pathway names. The gene (or gene-pair) symbol
# is only a few characters of a long shared template, so a fuzzy scorer rates
# wrong-gene siblings (e.g. "Signaling by BRCA1 mutants" vs "Signaling by AMER1
# mutants") highly. Each entry is (regex, capture-group-index); the captured
# token is the gene/key symbol used by the guard. Order matters — the first
# matching pattern wins, so the slash-pair and qualified variants precede the
# bare "signaling by <GENE>" catch-all.
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


def load_synonyms(reactome_ref: set) -> dict:
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


def fuzzy_canonical(
    pathway: str,
    reactome_ref: set,
    synonyms: dict,
    threshold: int = 85,
) -> tuple[bool, str, str]:
    """Resolve a pathway name to its canonical Reactome form.

    Returns ``(is_canonical, canonical_name, method)`` where ``method`` is one of
    ``exact`` / ``synonym`` / ``fuzzy`` / ``non_canonical``. On a canonical match
    ``canonical_name`` is the exact Reactome string; otherwise it is the bare
    input. Resolution proceeds in priority order:

      a) strip the R-HSA suffix, then a normalized exact match → "exact";
      b) the informal→canonical synonym map (case-insensitive) → "synonym";
      c) rapidfuzz token_sort_ratio ≥ ``threshold`` (gene-token guarded) → "fuzzy";
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
        return True, canonical, "exact"

    # b) Synonym map (lowercased keys; null entries already dropped on load).
    synonym = synonyms.get(bare.lower())
    if synonym:
        log.info("pathway %r → synonym %r", bare, synonym)
        return True, synonym, "synonym"

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
        if score >= threshold and gene_token_guard(bare, match):
            log.info("pathway %r → fuzzy %r (score=%.1f)", bare, match, score)
            return True, match, "fuzzy"
        if score >= threshold:
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
                threshold,
            )
    else:
        log.info("pathway %r → non_canonical (no reactome reference)", bare)

    # d) No canonical resolution.
    return False, bare, "non_canonical"
