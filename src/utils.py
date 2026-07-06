"""Shared utilities: logging, retry decorator, PMID validation, input loaders."""

from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path

import httpx
from tenacity import (
    before_sleep_log,
    retry as _tenacity_retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

LOG_DIR = Path("outputs")
LOG_FILE = LOG_DIR / "pipeline.log"

_LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"

PMID_RE = re.compile(r"^\d{7,8}$")

# Anthropic activates prompt caching only when the cached prefix is at least this
# many tokens (Opus/Sonnet). Used by the extractor/merger to sanity-check that
# their static system+tools prefix is large enough to cache.
CACHE_MIN_TOKENS = 1024


def estimate_tokens(text: str) -> int:
    """Rough token estimate (~4 chars/token), no network call.

    Used at import time to confirm a cacheable prefix clears CACHE_MIN_TOKENS
    without paying for an API count_tokens round-trip. JSON (tool schemas) tends
    to tokenize denser than this, so the estimate is conservative (under-counts),
    meaning a prefix that clears the threshold here clears it for real.
    """
    return len(text) // 4


def setup_logging(level: str = "INFO", log_dir: str | Path | None = None) -> logging.Logger:
    """Configure root logging to stdout and a ``pipeline.log`` file.

    ``log_dir`` chooses where the log file is written; it defaults to ``outputs/``
    for backward compatibility, but the pipelines pass their timestamped run
    directory so each run keeps its own log alongside its other artifacts.

    Returns the package logger ("bio_annot").
    """
    log_dir = Path(log_dir) if log_dir is not None else LOG_DIR
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "pipeline.log"

    formatter = logging.Formatter(_LOG_FORMAT)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)

    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Avoid duplicate handlers if called more than once.
    root.handlers.clear()
    root.addHandler(stream_handler)
    root.addHandler(file_handler)

    return logging.getLogger("bio_annot")


def resolve_run_dir() -> Path:
    """Directory of the run whose artifacts a *reader* tool should load.

    Runs write into a timestamped ``outputs/runs/<ts>/`` directory (see
    ``PipelineConfig.run_dir``) with ``outputs/latest`` pointing at the newest one.
    Reader tools (visualization, cytoscape export, synonym rebuild) resolve which
    run to read here, in priority order:

      1. ``$RUN_DIR`` if it names an existing directory — lets a caller pin a
         specific run (the pipeline sets this when it spawns build_synonyms.py so
         the child reads the exact run just written, not a fresh timestamp).
      2. ``outputs/latest`` — the most recent run.
      3. ``outputs/`` — legacy pre-timestamp layout / first-run fallback.
    """
    env_dir = os.getenv("RUN_DIR")
    if env_dir and Path(env_dir).is_dir():
        return Path(env_dir)
    latest = Path("outputs/latest")
    if latest.exists():
        return latest
    return Path("outputs")


# Module-level logger used by the retry decorator's before_sleep hook.
_retry_log = logging.getLogger("bio_annot.retry")

#: Retry decorator for HTTP-bound coroutines/functions.
#: 4 attempts, exponential backoff (2s..30s), retrying on httpx transport
#: and HTTP status errors, logging each retry.
retry = _tenacity_retry(
    stop=stop_after_attempt(4),
    wait=wait_exponential(min=2, max=30),
    retry=retry_if_exception_type((httpx.HTTPStatusError, httpx.RequestError)),
    before_sleep=before_sleep_log(_retry_log, logging.WARNING),
    reraise=True,
)


def validate_pmids(pmids: list[str]) -> list[str]:
    """Keep only strings that look like valid PMIDs (7-8 digits).

    Anything that does not match is logged and dropped. Never invents or
    passes through unvalidated PMIDs.
    """
    log = logging.getLogger("bio_annot.utils")
    valid: list[str] = []
    for pmid in pmids:
        s = str(pmid).strip()
        if PMID_RE.match(s):
            valid.append(s)
        else:
            log.warning("Dropping invalid PMID: %r", pmid)
    return valid


DEFAULT_DISEASE_CONTEXT = "cancer"
DEFAULT_DISEASE_TERMS = "cancer,tumor,carcinoma,sarcoma,lymphoma,leukemia"

# Appended to the disease-term-derived PubMed query as a generic catch-all so
# the search still surfaces disease-relevant abstracts beyond the listed terms.
_GENERIC_QUERY_TERM = "disease"
# How many of the (broad-first) DISEASE_TERMS to fold into the PubMed query,
# keeping the Entrez OR-clause focused rather than exhaustive.
_PUBMED_QUERY_TERM_COUNT = 2


def load_disease_context() -> dict:
    """Resolve the active disease context from the environment.

    Reads ``DISEASE_CONTEXT`` (a single label, e.g. "cancer") and
    ``DISEASE_TERMS`` (a comma-separated list, e.g.
    "cancer,tumor,carcinoma,..."), and returns:

      - ``context``: the DISEASE_CONTEXT label.
      - ``pubmed_query_terms``: a short, deduped, lowercase OR-clause for the
        PubMed search — the context plus the first couple of (broad-first)
        DISEASE_TERMS plus a generic "disease" catch-all. For the default
        cancer config this is ``["cancer", "tumor", "disease"]``.
      - ``scoring_terms``: the full DISEASE_TERMS set (plus the context),
        lowercased, used to match disease names during prioritization scoring.
    """
    context = (
        os.getenv("DISEASE_CONTEXT", DEFAULT_DISEASE_CONTEXT).strip()
        or DEFAULT_DISEASE_CONTEXT
    )
    raw_terms = os.getenv("DISEASE_TERMS", DEFAULT_DISEASE_TERMS)
    terms = [t.strip() for t in raw_terms.split(",") if t.strip()]
    if not terms:
        terms = [context]

    scoring_terms = {context.lower()} | {t.lower() for t in terms}

    # Build the focused PubMed OR-clause, deduped and order-preserving.
    pubmed_query_terms: list[str] = []
    seen: set[str] = set()
    for term in [context, *terms[:_PUBMED_QUERY_TERM_COUNT], _GENERIC_QUERY_TERM]:
        lowered = term.lower()
        if lowered not in seen:
            seen.add(lowered)
            pubmed_query_terms.append(lowered)

    return {
        "context": context,
        "pubmed_query_terms": pubmed_query_terms,
        "scoring_terms": scoring_terms,
    }


def load_gene_list(path: str) -> list[str]:
    """Read a gene-symbol file: strip whitespace, skip blanks and # comments.

    Returns uppercase gene symbols.
    """
    genes: list[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            genes.append(stripped.upper())
    return genes


def load_ref_set(path: str) -> set[str]:
    """Read a refs/*.txt file into a set of stripped, non-blank, non-comment lines."""
    refs: set[str] = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            refs.add(stripped)
    return refs
