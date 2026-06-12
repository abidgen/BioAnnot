"""STRING protein–protein interaction fetcher.

Pulls high-confidence interaction partners for a gene from the STRING database,
to enrich the target network with database-backed PPI edges that complement the
(sparser, less reliable) LLM-extracted interactors.
"""

from __future__ import annotations

import logging

import httpx

from src.utils import retry

log = logging.getLogger("bio_annot.string_db")

# interaction_partners is the correct endpoint for "partners of a single
# protein". The /network endpoint only returns edges *among* a set of input
# identifiers, so it yields nothing for a single gene.
STRING_PARTNERS_URL = "https://string-db.org/api/json/interaction_partners"

# STRING combined-score scale is 0–1000; 700 == "high confidence" per STRING docs.
DEFAULT_MIN_SCORE = 700
# Cap partners per gene so a hub protein doesn't flood the network.
DEFAULT_LIMIT = 50
# STRING asks API consumers to identify themselves on each request.
CALLER_IDENTITY = "bio_annotation_pipeline"


def _to_combined_score(raw) -> int | None:
    """Normalize STRING's score to the 0–1000 integer scale.

    STRING's JSON `score` is a 0–1 float (e.g. 0.999), while the API's
    required_score parameter is on the 0–1000 scale. Handle both conventions
    defensively: a value ≤ 1 is treated as a fraction and scaled up.
    """
    if raw is None:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return round(value * 1000) if value <= 1.0 else round(value)


@retry
async def fetch_string(
    gene: str,
    species: str = "9606",
    min_score: int = DEFAULT_MIN_SCORE,
    limit: int = DEFAULT_LIMIT,
) -> list[dict]:
    """Fetch high-confidence STRING interaction partners for a gene.

    Returns a list of ``{"partner": SYMBOL, "combined_score": int}`` (score on
    the 0–1000 scale), filtered to ``combined_score >= min_score``, deduped
    (keeping the max score per partner), and sorted by descending score.
    Partner symbols are uppercased to match the pipeline's gene symbols.
    Returns ``[]`` if the gene is unknown to STRING or has no qualifying partners.
    """
    params = {
        "identifiers": gene,
        "species": species,
        "required_score": str(min_score),
        "limit": str(limit),
        "caller_identity": CALLER_IDENTITY,
    }

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.get(STRING_PARTNERS_URL, params=params)
        # STRING returns 400 for an identifier it can't map; treat as "no partners".
        if resp.status_code == 400:
            log.warning("STRING: no mapping for gene %s (HTTP 400)", gene)
            return []
        resp.raise_for_status()
        rows = resp.json()

    # Keep the best score seen per partner (STRING may list a partner more than once).
    best: dict[str, int] = {}
    for row in rows:
        name = (row.get("preferredName_B") or "").upper()
        combined = _to_combined_score(row.get("score"))
        if not name or combined is None or combined < min_score:
            continue
        if name == gene.upper():
            continue
        best[name] = max(best.get(name, 0), combined)

    partners = [
        {"partner": name, "combined_score": score} for name, score in best.items()
    ]
    partners.sort(key=lambda p: p["combined_score"], reverse=True)

    log.info("STRING %s: %d partners (score >= %d)", gene, len(partners), min_score)
    return partners
