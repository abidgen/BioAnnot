"""UniProt REST fetcher: pull function, localization, keywords and GO terms for a gene."""

from __future__ import annotations

import logging
from typing import Any

import httpx

from src.utils import retry

log = logging.getLogger("bio_annot.uniprot")

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


@retry
async def fetch_uniprot(gene: str, organism: str = "9606") -> dict[str, Any]:
    """Fetch the top UniProtKB entry for a gene in a given organism (default human).

    Returns a flat dict with accession, names, function text, subcellular
    locations, keywords and GO terms. Returns {} if the gene is not found.
    """
    params = {
        "query": f"gene_exact:{gene} AND organism_id:{organism} AND reviewed:true",
        "format": "json",
        "size": "1",
    }

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.get(UNIPROT_SEARCH_URL, params=params)
        resp.raise_for_status()
        data = resp.json()

    results = data.get("results", [])
    if not results:
        log.warning("UniProt: no entry for gene %s (organism %s)", gene, organism)
        return {}

    entry = results[0]

    accession = entry.get("primaryAccession", "")
    gene_name = _extract_gene_name(entry)
    protein_name = _extract_protein_name(entry)
    function_text = _extract_function(entry)
    subcellular_locations = _extract_subcellular(entry)
    keywords = [kw.get("name", "") for kw in entry.get("keywords", []) if kw.get("name")]
    go_terms = _extract_go_terms(entry)

    result = {
        "accession": accession,
        "gene_name": gene_name,
        "protein_name": protein_name,
        "function_text": function_text,
        "subcellular_locations": subcellular_locations,
        "keywords": keywords,
        "go_terms": go_terms,
    }
    log.info(
        "UniProt %s: accession=%s, %d GO terms, %d locations",
        gene,
        accession,
        len(go_terms),
        len(subcellular_locations),
    )
    return result


def _extract_gene_name(entry: dict) -> str:
    genes = entry.get("genes", [])
    if genes:
        return genes[0].get("geneName", {}).get("value", "")
    return ""


def _extract_protein_name(entry: dict) -> str:
    return (
        entry.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")
    )


def _extract_function(entry: dict) -> str:
    """First FUNCTION comment's first text value."""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value", "")
    return ""


def _extract_subcellular(entry: dict) -> list[str]:
    """List of subcellular location value strings from SUBCELLULAR LOCATION comments."""
    locations: list[str] = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                value = loc.get("location", {}).get("value")
                if value:
                    locations.append(value)
    return locations


def _extract_go_terms(entry: dict) -> list[dict]:
    """GO cross-references → list of {id, aspect, term}.

    UniProt encodes GO aspect+term in a 'GoTerm' property like "C:nucleus",
    where the prefix (C/F/P) is the aspect and the remainder is the term.
    """
    go_terms: list[dict] = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        go_id = xref.get("id", "")
        aspect = ""
        term = ""
        for prop in xref.get("properties", []):
            if prop.get("key") == "GoTerm":
                raw = prop.get("value", "")
                if ":" in raw:
                    aspect, term = raw.split(":", 1)
                else:
                    term = raw
        go_terms.append({"id": go_id, "aspect": aspect, "term": term})
    return go_terms
