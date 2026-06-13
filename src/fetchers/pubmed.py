"""PubMed/Entrez fetcher: search PMIDs and fetch abstracts via NCBI E-utilities."""

from __future__ import annotations

import logging
import os

import httpx
import xmltodict

from src.utils import retry, validate_pmids, load_disease_context

log = logging.getLogger("bio_annot.pubmed")

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _entrez_params() -> dict:
    """Common Entrez params: email (NCBI policy) and optional API key."""
    params: dict[str, str] = {}
    email = os.getenv("NCBI_EMAIL")
    if email:
        params["email"] = email
    api_key = os.getenv("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    return params


def _as_list(value) -> list:
    """Normalize an xmltodict node that may be a dict, list, or None to a list."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


@retry
async def search_pmids(gene: str, max_results: int = 20) -> list[str]:
    """Search PubMed for PMIDs related to a gene in a disease context.

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
        **_entrez_params(),
    }

    async with httpx.AsyncClient(timeout=30.0) as client:
        resp = await client.get(ESEARCH_URL, params=params)
        resp.raise_for_status()
        data = resp.json()

    idlist = data.get("esearchresult", {}).get("idlist", [])
    pmids = validate_pmids(idlist)
    log.info("search_pmids(%s): %d PMIDs after validation", gene, len(pmids))
    return pmids


@retry
async def fetch_abstracts(pmids: list[str]) -> list[dict]:
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
