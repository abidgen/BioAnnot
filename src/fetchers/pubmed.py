"""PubMed/Entrez fetcher: search PMIDs and fetch abstracts via NCBI E-utilities."""

from __future__ import annotations

import logging
from typing import Any

import httpx
import xmltodict

from src.config import config
from src.utils import retry, validate_pmids, load_disease_context

log = logging.getLogger("bio_annot.pubmed")

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _entrez_params() -> dict:
    """Common Entrez params: email (NCBI policy) and optional API key."""
    params: dict[str, str] = {}
    if config.ncbi_email:
        params["email"] = config.ncbi_email
    if config.ncbi_api_key:
        params["api_key"] = config.ncbi_api_key
    return params


def _as_list(value) -> list:
    """Normalize an xmltodict node that may be a dict, list, or None to a list."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


@retry
async def search_pmids(
    gene: str,
    max_results: int = config.pubmed_max_results,
    limit: int = config.pubmed_extract_limit,
) -> list[str]:
    """Search PubMed for PMIDs related to a gene in a disease context.

    Fetches up to ``max_results`` candidates ranked by relevance (sort=relevance,
    not the ESearch default of most-recent), validates them, and returns the best
    ``limit``. The candidate pool is deliberately deeper than ``limit`` so that
    validate_pmids() drops don't erode the top relevance-ranked hits.

    Returns validated PMID strings (digits only, 7-8 chars).
    """
    query_terms = load_disease_context()["pubmed_query_terms"]
    or_clause = " OR ".join(query_terms)
    term = f"{gene}[gene] AND ({or_clause})"
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": str(max_results),
        "retmode": "json",
        # Rank by relevance so the abstracts that survive the extractor's
        # ~3000-word truncation are the most on-target, not merely the newest.
        "sort": "relevance",
        **_entrez_params(),
    }

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.get(ESEARCH_URL, params=params)
        resp.raise_for_status()
        data = resp.json()

    idlist = data.get("esearchresult", {}).get("idlist", [])
    pmids = validate_pmids(idlist)[:limit]
    log.info(
        "search_pmids(%s): %d PMIDs after validation (pool=%d, limit=%d)",
        gene, len(pmids), max_results, limit,
    )
    return pmids


@retry
async def fetch_abstracts(pmids: list[str]) -> list[dict[str, Any]]:
    """Fetch abstract records for a list of PMIDs.

    Returns a list of dicts: {pmid, title, abstract, year, journal}.
    Structured AbstractText (a list) is joined with spaces; missing abstracts
    yield an empty string and a warning.
    """
    if not pmids:
        return []

    params = {
        "db": "pubmed",
        "retmode": "xml",
        "id": ",".join(pmids),
        **_entrez_params(),
    }

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.get(EFETCH_URL, params=params)
        resp.raise_for_status()
        parsed = xmltodict.parse(resp.text)

    articles = _as_list(
        parsed.get("PubmedArticleSet", {}).get("PubmedArticle")
    )

    results: list[dict] = []
    for art in articles:
        citation = art.get("MedlineCitation", {})
        pmid = _extract_pmid(citation)
        article = citation.get("Article", {})

        title = _extract_text(article.get("ArticleTitle"))
        abstract = _extract_abstract(article)
        year = _extract_year(article)
        journal = _extract_text(article.get("Journal", {}).get("Title"))

        if not abstract:
            log.warning("No abstract for PMID %s", pmid)

        results.append(
            {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "year": year,
                "journal": journal,
            }
        )

    log.info("fetch_abstracts: parsed %d records", len(results))
    return results


def _extract_pmid(citation: dict) -> str:
    """PMID may be a plain string or a dict with a '#text' key."""
    pmid = citation.get("PMID")
    if isinstance(pmid, dict):
        return str(pmid.get("#text", ""))
    return str(pmid) if pmid is not None else ""


def _extract_text(node) -> str:
    """Extract text from a node that may be a string or a dict with '#text'."""
    if node is None:
        return ""
    if isinstance(node, dict):
        return str(node.get("#text", "")).strip()
    return str(node).strip()


def _extract_abstract(article: dict) -> str:
    """Join AbstractText nodes (handles structured abstracts as a list)."""
    abstract_node = article.get("Abstract", {})
    if not abstract_node:
        return ""
    texts = _as_list(abstract_node.get("AbstractText"))
    parts = [_extract_text(t) for t in texts]
    return " ".join(p for p in parts if p)


def _extract_year(article: dict) -> str:
    """Pull publication year from Journal > JournalIssue > PubDate."""
    pub_date = (
        article.get("Journal", {})
        .get("JournalIssue", {})
        .get("PubDate", {})
    )
    if not isinstance(pub_date, dict):
        return ""
    if pub_date.get("Year"):
        return str(pub_date["Year"])
    # MedlineDate is a free-text fallback, e.g. "2019 Jan-Feb".
    medline = pub_date.get("MedlineDate", "")
    return str(medline)[:4] if medline else ""
