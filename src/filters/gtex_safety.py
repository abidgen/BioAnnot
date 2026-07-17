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
from typing import Any

import httpx
import pandas as pd

from src.config import config

log = logging.getLogger("bio_annot.gtex_safety")

# GTEx moved its v8 bulk downloads to the "adult-gtex" bucket; the older
# gtex_analysis_v8/rna_seq_data/ path now 404s. Same filename, new location.
GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)
CACHE_PATH = Path("refs/gtex_median_tpm.gct.gz")

# Two-tier normal-tissue safety model. Names must match GTEx column headers exactly.
#
# TIER 1 — vital organs: expression in any ONE of these above the (low) vital
# threshold is a hard safety flag; damaging these tissues is poorly tolerated.
TIER1_VITAL_TISSUES = {
    "Brain - Cortex",
    "Brain - Cerebellum",
    "Heart - Left Ventricle",
    "Liver",
    "Kidney - Cortex",
    "Lung",
    "Adrenal Gland",
}

# TIER 2 — secondary sensitive tissues: a soft flag requires expression in at
# least ``gtex_tier2_min_tissues`` of these above the (higher) standard threshold.
TIER2_SENSITIVE_TISSUES = {
    "Small Intestine - Terminal Ileum",
    "Colon - Sigmoid",
    "Spleen",
    "Skin - Sun Exposed (Lower leg)",
    "Whole Blood",
}

# All tissues considered by the safety filter (union of both tiers).
SENSITIVE_TISSUES = TIER1_VITAL_TISSUES | TIER2_SENSITIVE_TISSUES

# Neutral multiplier for an unflagged target. The tier-1 and tier-2 penalties are
# config-driven (config.gtex_tier1_penalty / gtex_tier2_penalty); tier 1 takes
# precedence over tier 2.
NO_PENALTY = 1.00

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
        "tier1_flag": False,
        "tier2_flag": False,
        "tier1_high_tissues": {},
        "tier2_high_tissues": {},
        "high_expression_tissues": [],
        "max_vital_tpm": 0.0,
        "max_tpm": 0.0,
        "tissue_count_above_threshold": 0,
        "safety_penalty": NO_PENALTY,
        "safety_note": "",
    }


def assess_safety(
    gene: str,
    vital_tpm_threshold: float = config.gtex_vital_tpm_threshold,
    tpm_threshold: float = config.gtex_tpm_threshold,
    tier2_min_tissues: int = config.gtex_tier2_min_tissues,
    tier1_penalty: float = config.gtex_tier1_penalty,
    tier2_penalty: float = config.gtex_tier2_penalty,
) -> dict[str, Any]:
    """Assess on-target normal-tissue safety risk for a gene from GTEx (two-tier).

    TIER 1 (vital organs): expression above ``vital_tpm_threshold`` in ANY single
    tier-1 tissue is a hard flag (``tier1_flag``), penalty ``tier1_penalty``.
    TIER 2 (secondary sensitive tissues): expression above ``tpm_threshold`` in at
    least ``tier2_min_tissues`` tier-2 tissues is a soft flag (``tier2_flag``),
    penalty ``tier2_penalty``. Tier 1 takes precedence when both fire.

    A flag deprioritizes the target downstream (score × penalty), never eliminates
    it. Returns an unflagged, no-penalty assessment if the gene is not in GTEx.
    """
    if _GTEX is None or gene not in _GTEX.index:
        return _empty_assessment(gene)

    row = _GTEX.loc[gene]
    # A symbol can map to multiple rows (paralogous Ensembl IDs); collapse by the
    # max median TPM per tissue so we assess the worst case.
    if isinstance(row, pd.DataFrame):
        row = row.max(numeric_only=True)

    # Tier 1 — vital organs (hard flag on any single organ over the vital threshold).
    tier1_tpm = {t: float(row[t]) for t in TIER1_VITAL_TISSUES if t in row.index}
    tier1_high_tissues = {
        t: tpm for t, tpm in tier1_tpm.items() if tpm > vital_tpm_threshold
    }
    tier1_flag = len(tier1_high_tissues) >= 1
    max_vital_tpm = max(tier1_tpm.values()) if tier1_tpm else 0.0

    # Tier 2 — secondary tissues (soft flag once enough clear the standard threshold).
    tier2_tpm = {t: float(row[t]) for t in TIER2_SENSITIVE_TISSUES if t in row.index}
    tier2_high_tissues = {
        t: tpm for t, tpm in tier2_tpm.items() if tpm > tpm_threshold
    }
    tier2_flag = len(tier2_high_tissues) >= tier2_min_tissues

    # Tier 1 penalty takes precedence over tier 2.
    if tier1_flag:
        safety_penalty = tier1_penalty
    elif tier2_flag:
        safety_penalty = tier2_penalty
    else:
        safety_penalty = NO_PENALTY

    safety_flag = tier1_flag or tier2_flag
    # Backward-compatible union field: all tissues over their tier's threshold.
    high_expression_tissues = sorted(
        set(tier1_high_tissues) | set(tier2_high_tissues)
    )
    all_tpm = {**tier1_tpm, **tier2_tpm}
    max_tpm = max(all_tpm.values()) if all_tpm else 0.0

    if tier1_flag:
        note = (
            f"Vital-organ expression over {vital_tpm_threshold} TPM in "
            f"{sorted(tier1_high_tissues)} — hard safety flag. {TF_CONTEXT_CAVEAT}"
        )
    elif tier2_flag:
        note = (
            f"Elevated expression in {len(tier2_high_tissues)} secondary sensitive "
            f"tissues — soft safety flag. {TF_CONTEXT_CAVEAT}"
        )
    else:
        note = ""

    return {
        "gene": gene,
        "safety_flag": safety_flag,
        "tier1_flag": tier1_flag,
        "tier2_flag": tier2_flag,
        "tier1_high_tissues": tier1_high_tissues,
        "tier2_high_tissues": tier2_high_tissues,
        "high_expression_tissues": high_expression_tissues,
        "max_vital_tpm": max_vital_tpm,
        "max_tpm": max_tpm,
        "tissue_count_above_threshold": len(high_expression_tissues),
        "safety_penalty": safety_penalty,
        "safety_note": note,
    }
