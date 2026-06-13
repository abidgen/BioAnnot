"""CellxGene Census single-cell expression fetcher.

Grounds the LLM-extracted ``cellular_states`` field in measured per-cell-type
expression: for a gene and a broad tissue (default lung), it returns the mean
expression in each cell type, restricted to cell types with enough cells to be
trustworthy.

Heavy lifting (opening the S3-backed census, scanning cell metadata) runs in a
worker thread via ``asyncio.to_thread`` so the async pipeline's event loop is
not blocked. The expensive query result is cached to disk per (gene, tissue);
the ``cellxgene_census`` import is deferred into the query path so importing this
module never requires the package to be installed (e.g. when ENABLE_CELLXGENE
is off).
"""

from __future__ import annotations

import asyncio
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from src.config import config
from tenacity import (
    before_sleep_log,
    retry as _tenacity_retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

log = logging.getLogger("bio_annot.cellxgene")

# The shared utils.retry only catches httpx errors, but the census stack fails
# through TileDB/S3 (botocore, tiledbsoma), not httpx. Use a census-specific
# retry that catches any exception, backs off longer (S3 throttling), and logs
# each attempt with the exception type via before_sleep_log.
_retry_log = logging.getLogger("bio_annot.cellxgene.retry")
census_retry = _tenacity_retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(min=5, max=60),
    retry=retry_if_exception_type(Exception),
    before_sleep=before_sleep_log(_retry_log, logging.WARNING),
    reraise=True,
)

# Pin the census version for reproducibility — "latest" drifts as new releases
# land, which would silently change results between runs. Centralized in
# src.config (env-configurable).
CENSUS_VERSION = config.census_version
CENSUS_TISSUE = config.census_tissue
CENSUS_MIN_CELLS = config.census_min_cells
CACHE_DIR = Path(config.census_cache_dir)
ENABLE_CELLXGENE = config.enable_cellxgene

CENSUS_ORGANISM = "Homo sapiens"


def _cache_path(gene: str, tissue: str) -> Path:
    """Per-(gene, tissue) cache file. min_cells is applied after load, not here,
    so changing CENSUS_MIN_CELLS reuses one cache file rather than refetching."""
    return CACHE_DIR / f"{gene}_{tissue}.json"


def _load_cache(gene: str, tissue: str) -> dict | None:
    """Return cached per-cell-type stats, or None if absent/unreadable."""
    path = _cache_path(gene, tissue)
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("CellxGene cache unreadable (%s); refetching: %s", exc, path)
        return None


def _write_cache(gene: str, tissue: str, stats: dict) -> None:
    """Atomically write per-cell-type stats to the cache."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = _cache_path(gene, tissue)
    tmp = path.with_suffix(".json.part")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)
    tmp.rename(path)  # atomic: a partial write never looks complete


def _get_anndata(census, gene: str, tissue: str):
    """Slice the census to one gene in one tissue (primary cells only).

    get_anndata's column-selection kwarg was renamed across cellxgene_census
    releases (``column_names`` dict → ``obs_column_names``/``var_column_names``),
    so try the newer form and fall back to the older one.
    """
    import cellxgene_census

    common = dict(
        organism=CENSUS_ORGANISM,
        measurement_name="RNA",
        X_name="raw",
        var_value_filter=f"feature_name == '{gene}'",
        # is_primary_data == True drops cells duplicated across datasets, so a
        # cell is never counted twice in the per-cell-type mean.
        obs_value_filter=(
            f"tissue_general == '{tissue}' and is_primary_data == True"
        ),
    )
    try:
        return cellxgene_census.get_anndata(
            census,
            **common,
            obs_column_names=["cell_type"],
            var_column_names=["feature_name"],
        )
    except TypeError:
        return cellxgene_census.get_anndata(
            census,
            **common,
            column_names={"obs": ["cell_type"], "var": ["feature_name"]},
        )


@census_retry
def _query_census(gene: str, tissue: str) -> dict:
    """Blocking census query → ``{cell_type: {"mean": float, "n_cells": int}}``.

    Computes the per-cell-type mean of raw counts and the cell count, unfiltered
    (the caller applies min_cells). Returns an empty dict if the gene/tissue
    combination has no cells. Runs in a worker thread — never call directly from
    the event loop.
    """
    import cellxgene_census  # noqa: F401 — defer import; presence checked here

    log.info(
        "CellxGene: querying census %s for %s in tissue_general=%s",
        CENSUS_VERSION,
        gene,
        tissue,
    )
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        adata = _get_anndata(census, gene, tissue)

    if adata.n_obs == 0 or adata.n_vars == 0:
        log.warning("CellxGene: no cells for %s in %s", gene, tissue)
        return {}

    # adata.X is cells × 1 (single gene); densify the column and group by type.
    matrix = adata.X
    expr = (
        np.asarray(matrix.todense()).ravel()
        if hasattr(matrix, "todense")
        else np.asarray(matrix).ravel()
    )
    frame = pd.DataFrame(
        {"cell_type": adata.obs["cell_type"].to_numpy(), "expr": expr}
    )
    agg = frame.groupby("cell_type", observed=True)["expr"].agg(["mean", "size"])

    stats = {
        str(cell_type): {"mean": float(row["mean"]), "n_cells": int(row["size"])}
        for cell_type, row in agg.iterrows()
    }
    log.info(
        "CellxGene: %s/%s → %d cell types, %d cells total",
        gene,
        tissue,
        len(stats),
        int(adata.n_obs),
    )
    return stats


async def fetch_cellxgene(
    gene: str,
    tissue: str | None = None,
    min_cells: int | None = None,
) -> dict[str, float]:
    """Mean per-cell-type expression for a gene in a tissue, from CellxGene Census.

    Returns ``{cell_type: mean_expr}`` (raw-count mean), restricted to cell types
    with at least ``min_cells`` cells and sorted by descending mean expression.
    Results are cached to ``CENSUS_CACHE_DIR/{gene}_{tissue}.json``. Returns ``{}``
    if the fetcher is disabled (ENABLE_CELLXGENE) or the gene/tissue has no data.
    """
    tissue = tissue or CENSUS_TISSUE
    min_cells = CENSUS_MIN_CELLS if min_cells is None else min_cells
    gene = gene.upper()  # census feature_name uses HGNC symbols

    if not ENABLE_CELLXGENE:
        log.info("CellxGene disabled (ENABLE_CELLXGENE); skipping %s", gene)
        return {}

    stats = _load_cache(gene, tissue)
    if stats is None:
        stats = await asyncio.to_thread(_query_census, gene, tissue)
        _write_cache(gene, tissue, stats)

    filtered = {
        cell_type: s["mean"]
        for cell_type, s in stats.items()
        if s.get("n_cells", 0) >= min_cells
    }
    return dict(
        sorted(filtered.items(), key=lambda kv: kv[1], reverse=True)
    )
