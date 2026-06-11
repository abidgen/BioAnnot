"""Shared utilities: logging, retry decorator, PMID validation, input loaders."""

from __future__ import annotations

import logging
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


def setup_logging(level: str = "INFO") -> logging.Logger:
    """Configure root logging to stdout and outputs/pipeline.log.

    Returns the package logger ("bio_annot").
    """
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    formatter = logging.Formatter(_LOG_FORMAT)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(LOG_FILE)
    file_handler.setFormatter(formatter)

    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Avoid duplicate handlers if called more than once.
    root.handlers.clear()
    root.addHandler(stream_handler)
    root.addHandler(file_handler)

    return logging.getLogger("bio_annot")


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
