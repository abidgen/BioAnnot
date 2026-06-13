"""Build/update the informal→canonical pathway synonym map (local-first).

For every NON-CANONICAL pathway name emitted into outputs/final_annotations.json,
resolve a canonical Reactome name with a cheap local-first cascade, only falling
back to an LLM for the genuinely ambiguous names:

  1. exact  — normalized exact match against the full refs/reactome_pathways.txt
              (no API call)
  2. fuzzy  — rapidfuzz token_sort_ratio over the FULL reference set; accept the
              best hit if its score >= FUZZY_THRESHOLD (no API call)
  3. llm    — names that fail steps 1 and 2 go to SYNONYM_MODEL in one call,
              each with its top-SYNONYM_CANDIDATES rapidfuzz candidates; the model
              picks one candidate or null, and the choice is validated against the
              reference before being kept (method "llm" or "null")

The result is merged into refs/pathway_synonyms.json and consumed by
src.merger's fuzzy canonicalization. Incremental: names already in the map are
skipped. Run standalone, or via AUTO_UPDATE_SYNONYMS after a pipeline run:

    python scripts/build_synonyms.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

from rapidfuzz import fuzz, process

# Make the project root importable when run as `python scripts/build_synonyms.py`
# (the script's own dir, not the project root, is on sys.path by default).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.config import config
from src.pathways import (
    NON_CANONICAL_PREFIX,
    _normalize,
    extract_key_token,
    gene_token_guard,
)

FINAL_ANNOTATIONS = Path("outputs/final_annotations.json")
SYNONYMS_PATH = Path("refs/pathway_synonyms.json")
REACTOME_PATH = Path("refs/reactome_pathways.txt")

# Model used only for step 3; both thresholds are centralized in src.config.
SYNONYM_MODEL = config.synonym_model
FUZZY_THRESHOLD = config.fuzzy_threshold
SYNONYM_CANDIDATES = config.synonym_candidates

# Pathway normalization, the templated-pattern list, and the gene-token guard are
# imported from src.pathways (shared with src.merger), so there is one copy.


def load_noncanonical_names() -> list[str]:
    """Collect distinct bare NON-CANONICAL pathway names from the annotations."""
    if not FINAL_ANNOTATIONS.exists():
        print(f"No {FINAL_ANNOTATIONS} — run the pipeline first.")
        return []
    with open(FINAL_ANNOTATIONS, encoding="utf-8") as f:
        data = json.load(f)
    names: set[str] = set()
    for record in data.values():
        for pathway in record.get("pathways", []):
            if pathway.startswith(NON_CANONICAL_PREFIX):
                names.add(pathway[len(NON_CANONICAL_PREFIX):].strip())
    return sorted(names)


def load_synonyms() -> dict:
    """Load the existing synonym map, or an empty map if absent/unreadable."""
    if not SYNONYMS_PATH.exists():
        return {}
    try:
        with open(SYNONYMS_PATH, encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        print(f"Could not read {SYNONYMS_PATH} ({exc}); starting fresh.")
        return {}


def load_reactome_ref() -> list[str]:
    """Load the full canonical Reactome name list (empty if the file is missing)."""
    if not REACTOME_PATH.exists():
        return []
    return [
        line.strip()
        for line in REACTOME_PATH.read_text().splitlines()
        if line.strip()
    ]


def _parse_json_object(text: str) -> dict:
    """Extract a JSON object from a model response (tolerates code fences)."""
    text = text.strip()
    if text.startswith("```"):
        text = text.split("```", 2)[1]
        if text.lower().startswith("json"):
            text = text[4:]
    start, end = text.find("{"), text.rfind("}")
    if start == -1 or end == -1:
        raise ValueError(f"No JSON object found in model response: {text[:200]!r}")
    return json.loads(text[start : end + 1])


def query_llm(candidates_by_name: dict[str, list[str]]) -> dict:
    """Ask the model to pick the best canonical name from each name's candidates."""
    import anthropic

    client = anthropic.Anthropic()
    system = (
        "You are a pathway nomenclature expert. For each informal pathway name you "
        "are given a list of candidate canonical Reactome pathway names. Choose the "
        "single candidate that is the true canonical equivalent, or null if none of "
        "the candidates is a genuine match. Choose only from the provided "
        "candidates. Output strictly valid JSON and nothing else."
    )
    user = (
        "Map each informal pathway name to the best canonical Reactome name from "
        "its candidate list (or null if none is a true match):\n\n"
        f"{json.dumps(candidates_by_name, indent=2)}\n\n"
        "Return JSON only: {informal_name: chosen_candidate_or_null}"
    )
    response = client.messages.create(
        model=SYNONYM_MODEL,
        max_tokens=4096,
        system=system,
        messages=[{"role": "user", "content": user}],
    )
    usage = response.usage
    print(
        f"LLM call ({SYNONYM_MODEL}): input_tokens={usage.input_tokens} "
        f"output_tokens={usage.output_tokens}"
    )
    text = "".join(
        b.text for b in response.content if getattr(b, "type", None) == "text"
    )
    return _parse_json_object(text)


def main() -> None:
    noncanonical = load_noncanonical_names()
    if not noncanonical:
        print("No NON-CANONICAL pathway names found; nothing to do.")
        return

    synonyms = load_synonyms()
    new_names = [name for name in noncanonical if name not in synonyms]
    print(
        f"Found {len(noncanonical)} distinct NON-CANONICAL names "
        f"({len(synonyms)} already mapped, {len(new_names)} new)."
    )
    if not new_names:
        print("Synonym map already covers every NON-CANONICAL name; up to date.")
        return

    reactome_ref = load_reactome_ref()
    canonical_by_norm = {_normalize(r): r for r in reactome_ref}
    print(f"Loaded {len(reactome_ref)} canonical Reactome names.\n")

    # results[name] = (canonical_or_None, method, score_or_None)
    results: dict[str, tuple] = {}
    candidates_by_name: dict[str, list[str]] = {}
    # rejected[name] = (candidate, q_token, c_token, score) for a top fuzzy hit
    # that cleared the score but was vetoed by the gene-token guard.
    rejected: dict[str, tuple] = {}

    for name in new_names:
        # Step 1 — exact normalized match (no API).
        canonical = canonical_by_norm.get(_normalize(name))
        if canonical is not None:
            results[name] = (canonical, "exact", 100.0)
            continue

        # Step 2 — fuzzy over the FULL reference (no API). Pull a wide candidate
        # list (3× the LLM budget) so guard-vetoed siblings can be skipped while
        # still leaving enough genuine options for the LLM step.
        ranked = process.extract(
            name,
            reactome_ref,
            scorer=fuzz.token_sort_ratio,
            processor=str.lower,
            limit=SYNONYM_CANDIDATES * 3,
        )

        # Accept the highest-scoring candidate that clears the threshold AND
        # passes the gene-token guard. A guard-vetoed high scorer is logged and
        # skipped (demoted toward the LLM) rather than accepted.
        fuzzy_match = None
        for cand, score, _ in ranked:
            if score < FUZZY_THRESHOLD:
                break
            if gene_token_guard(name, cand):
                fuzzy_match = (cand, score)
                break
            q_token, c_token = extract_key_token(name), extract_key_token(cand)
            print(
                f"  Gene token mismatch: {name} → {cand} "
                f"({q_token} ≠ {c_token}), demoting to LLM"
            )
            rejected.setdefault(name, (cand, q_token, c_token, score))

        if fuzzy_match is not None:
            results[name] = (fuzzy_match[0], "fuzzy", fuzzy_match[1])
            continue

        # Step 3 — defer to the LLM. Candidates are guard-filtered (wrong-gene
        # siblings removed) so the model only ever sees plausible options.
        candidates_by_name[name] = [
            cand for cand, _score, _idx in ranked if gene_token_guard(name, cand)
        ][:SYNONYM_CANDIDATES]
        results[name] = (None, "pending", ranked[0][1] if ranked else 0.0)

    # One LLM call for everything that failed local resolution.
    if candidates_by_name:
        print(
            f"{len(candidates_by_name)} name(s) need the LLM "
            f"(top-{SYNONYM_CANDIDATES} candidates each)…"
        )
        mappings = query_llm(candidates_by_name)
        for name in candidates_by_name:
            chosen = mappings.get(name)
            # Validate the model's choice against the reference before trusting it.
            canonical = canonical_by_norm.get(_normalize(chosen)) if chosen else None
            results[name] = (
                (canonical, "llm", None) if canonical else (None, "null", None)
            )

    # Merge into the map (canonical name, or null for unresolved) and save.
    for name, (canonical, _method, _score) in results.items():
        synonyms[name] = canonical
    SYNONYMS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(SYNONYMS_PATH, "w", encoding="utf-8") as f:
        json.dump(synonyms, f, indent=2, sort_keys=True)

    by_method = Counter(method for _c, method, _s in results.values())
    print(
        f"\nResolution: exact={by_method['exact']} fuzzy={by_method['fuzzy']} "
        f"llm={by_method['llm']} null={by_method['null']} "
        f"(LLM calls: {1 if candidates_by_name else 0})"
    )
    print("Per-name resolution:")
    for name in new_names:
        canonical, method, score = results[name]
        if method == "exact":
            print(f"  {name!r} → exact → {canonical!r}")
            continue
        if method == "fuzzy":
            print(f"  {name!r} → fuzzy (score={score:.0f}) → {canonical!r}")
            continue
        # llm / null — show the guard rejection that demoted it, if any.
        rej = ""
        if name in rejected:
            _cand, q_token, c_token, rscore = rejected[name]
            rej = f"fuzzy REJECTED ({q_token} ≠ {c_token}, score={rscore:.0f}) → "
        target = repr(canonical) if method == "llm" else "null"
        print(f"  {name!r} → {rej}LLM → {target}")
    print(f"\nWrote {len(synonyms)} total entries → {SYNONYMS_PATH}")


if __name__ == "__main__":
    main()
