"""OpenTargets GraphQL fetcher: resolve a gene symbol and pull associations."""

from __future__ import annotations

import logging

import httpx

from src.utils import retry

log = logging.getLogger("bio_annot.opentargets")

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

_SEARCH_QUERY = """
query ($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id }
  }
}
"""

_TARGET_QUERY = """
query ($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    functionDescriptions
    pathways { pathway pathwayId }
    associatedDiseases(page: {index: 0, size: 10}) {
      rows {
        disease { name id }
        score
      }
    }
  }
}
"""


async def _graphql(client: httpx.AsyncClient, query: str, variables: dict) -> dict:
    """POST a GraphQL query and return the 'data' object."""
    resp = await client.post(
        OPENTARGETS_GRAPHQL,
        json={"query": query, "variables": variables},
    )
    resp.raise_for_status()
    payload = resp.json()
    if payload.get("errors"):
        log.warning("OpenTargets GraphQL errors: %s", payload["errors"])
    return payload.get("data") or {}


@retry
async def fetch_opentargets(gene_symbol: str) -> dict:
    """Resolve a gene symbol to an Ensembl target ID, then fetch associations.

    Returns a dict with ensembl_id, symbol, function_descriptions, pathways and
    disease_associations. Returns {} if the symbol cannot be resolved.
    """
    async with httpx.AsyncClient(timeout=30.0) as client:
        # 1. Resolve gene symbol -> Ensembl ID.
        search_data = await _graphql(client, _SEARCH_QUERY, {"q": gene_symbol})
        hits = search_data.get("search", {}).get("hits", [])
        if not hits:
            log.warning("OpenTargets: no target hit for %s", gene_symbol)
            return {}
        ensembl_id = hits[0].get("id", "")
        if not ensembl_id:
            log.warning("OpenTargets: empty Ensembl ID for %s", gene_symbol)
            return {}

        # 2. Fetch target associations.
        target_data = await _graphql(
            client, _TARGET_QUERY, {"ensemblId": ensembl_id}
        )

    target = target_data.get("target")
    if not target:
        log.warning("OpenTargets: no target record for %s (%s)", gene_symbol, ensembl_id)
        return {}

    pathways = [
        {"pathway": p.get("pathway", ""), "pathwayId": p.get("pathwayId", "")}
        for p in (target.get("pathways") or [])
    ]

    disease_rows = (target.get("associatedDiseases") or {}).get("rows") or []
    disease_associations = [
        {
            "name": row.get("disease", {}).get("name", ""),
            "id": row.get("disease", {}).get("id", ""),
            "score": row.get("score"),
        }
        for row in disease_rows
    ]

    result = {
        "ensembl_id": target.get("id", ensembl_id),
        "symbol": target.get("approvedSymbol", gene_symbol),
        "approved_name": target.get("approvedName", ""),
        "function_descriptions": target.get("functionDescriptions") or [],
        "pathways": pathways,
        "disease_associations": disease_associations,
    }
    log.info(
        "OpenTargets %s: %s, %d pathways, %d diseases",
        gene_symbol,
        result["ensembl_id"],
        len(pathways),
        len(disease_associations),
    )
    return result
