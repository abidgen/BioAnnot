"""GTEx normal-tissue safety filter.

Flags candidate targets that are highly expressed in sensitive normal tissues
(a classic on-target/off-tumor safety concern). Expression comes from the GTEx
v8 gene-level median-TPM table, downloaded and cached once on first import.

The downloaded table is large (~7 MB gz, ~56k genes × 54 tissues), so it is
parsed once at module load into a module-level DataFrame indexed by gene symbol.
"""

from __future__ import annotations

import logging
from pathlib import Path

import httpx
import pandas as pd

log = logging.getLogger("bio_annot.gtex_safety")

# GTEx moved its v8 bulk downloads to the "adult-gtex" bucket; the older
# gtex_analysis_v8/rna_seq_data/ path now 404s. Same filename, new location.
GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)
CACHE_PATH = Path("refs/gtex_median_tpm.gct.gz")

# Sensitive normal tissues: high expression here raises a safety concern for a
# therapeutic target. Names must match GTEx column headers exactly.
SENSITIVE_TISSUES = {
    "Brain - Cortex",
    "Brain - Cerebellum",
    "Heart - Left Ventricle",
    "Liver",
    "Kidney - Cortex",
    "Lung",
    "Adrenal Gland",
    "Small Intestine - Terminal Ileum",
    "Colon - Sigmoid",
}

# Appended to the safety_note of flagged genes: a high-expression flag is a
# prompt for review, not a verdict — for tissue-specific factors it can reflect
# on-target biology rather than a true off-tumor liability.
TF_CONTEXT_CAVEAT = (
    "For tissue-specific transcription factors (e.g. FOXF1 in lung), high normal "
    "tissue expression may reflect on-target biology rather than a true safety "
    "liability. Manual review recommended for safety_flag=True genes."
)


def _download_gtex(url: str, dest: Path) -> None:
    """Stream the GTEx median-TPM table to ``dest`` (idempotent: skips if present)."""
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    log.info("Downloading GTEx median-TPM table → %s", dest)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with httpx.stream("GET", url, timeout=120.0, follow_redirects=True) as resp:
        resp.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in resp.iter_bytes():
                f.write(chunk)
    tmp.rename(dest)  # atomic: a partial download never looks complete


def _load_gtex() -> pd.DataFrame | None:
    """Download (if needed) and parse the GTEx GCT into a symbol-indexed frame.

    GCT format: line 1 is "#1.2", line 2 is "<nrows>\\t<ncols>", line 3 is the
    header (Name, Description, then one column per tissue). Returns None if the
    table can't be obtained or parsed, so importing this module never breaks the
    pipeline when offline — assess_safety then degrades to "no concern".
    """
    try:
        _download_gtex(GTEX_URL, CACHE_PATH)
        frame = pd.read_csv(CACHE_PATH, sep="\t", skiprows=2, compression="gzip")
        frame = frame.set_index("Description")
    except Exception as exc:  # noqa: BLE001 — never let a load failure break import
        log.error("Could not load GTEx table (%s); safety filter disabled", exc)
        return None

    missing = sorted(SENSITIVE_TISSUES - set(frame.columns))
    if missing:
        log.warning("GTEx table missing expected sensitive-tissue columns: %s", missing)
    log.info(
        "Loaded GTEx median-TPM table: %d genes × %d tissues",
        frame.shape[0],
        frame.shape[1],
    )
    return frame


# Module-level cache: parsed once on first import.
_GTEX: pd.DataFrame | None = _load_gtex()


def _empty_assessment(gene: str) -> dict:
    """Safe default when a gene is absent from GTEx (or the table is unavailable)."""
    return {
        "gene": gene,
        "safety_flag": False,
        "high_expression_tissues": [],
        "max_tpm": 0.0,
        "tissue_count_above_threshold": 0,
        "safety_note": "",
    }


def assess_safety(
    gene: str, tpm_threshold: float = 10.0, min_tissues: int = 3
) -> dict:
    """Assess on-target normal-tissue safety risk for a gene from GTEx.

    Counts sensitive tissues whose median TPM exceeds ``tpm_threshold``; if at
    least ``min_tissues`` qualify, the gene is flagged as a potential safety
    concern (deprioritized downstream, not eliminated). ``max_tpm`` is the
    highest median TPM across the sensitive tissues. Returns an unflagged
    assessment with empty lists if the gene is not in GTEx.
    """
    if _GTEX is None or gene not in _GTEX.index:
        return _empty_assessment(gene)

    row = _GTEX.loc[gene]
    # A symbol can map to multiple rows (paralogous Ensembl IDs); collapse by the
    # max median TPM per tissue so we assess the worst case.
    if isinstance(row, pd.DataFrame):
        row = row.max(numeric_only=True)

    present_tissues = [t for t in SENSITIVE_TISSUES if t in row.index]
    tissue_tpm = {t: float(row[t]) for t in present_tissues}

    high_expression_tissues = sorted(
        t for t, tpm in tissue_tpm.items() if tpm > tpm_threshold
    )
    count = len(high_expression_tissues)
    max_tpm = max(tissue_tpm.values()) if tissue_tpm else 0.0
    flagged = count >= min_tissues

    note = (
        f"High normal tissue expression in {count} sensitive tissues — "
        f"review before advancing. {TF_CONTEXT_CAVEAT}"
        if flagged
        else ""
    )

    return {
        "gene": gene,
        "safety_flag": flagged,
        "high_expression_tissues": high_expression_tissues,
        "max_tpm": max_tpm,
        "tissue_count_above_threshold": count,
        "safety_note": note,
    }
