"""Tests for run_gene's two-layer cache execution order and stat accounting.

run_gene chains fetch → extract → merge → enrich behind a two-layer cache:
  1. final-cache hit  → return the record, run nothing (genes_cached++)
  2. raw-cache hit     → skip fetch + extract, replay merge + enrich (genes_remerged++)
  3. miss on both      → run the whole chain, write the raw cache

The stage functions and cache I/O are replaced with mocks, so no network, API, or
real disk is touched (only the annotations log, which is redirected to tmp_path).

Importing ``pipeline`` runs ``load_dotenv()`` at import; snapshot/restore
``os.environ`` so it has no global side effect (keeps test_config.py's clean-env
assertions valid).
"""

import os as _os
from unittest.mock import AsyncMock, Mock

import pytest

_ENV_SNAPSHOT = dict(_os.environ)
import pipeline  # noqa: E402
from pipeline import run_gene, RunStats  # noqa: E402
_os.environ.clear()
_os.environ.update(_ENV_SNAPSHOT)

from src.config import PipelineConfig  # noqa: E402


@pytest.fixture
def config():
    return PipelineConfig()


@pytest.fixture
def stats():
    return RunStats()


def _patch_stages(monkeypatch):
    """Replace the four stage functions with AsyncMocks; return them."""
    fetch = AsyncMock(return_value={
        "pubmed_text": "", "pmids": [], "uniprot": {}, "opentargets": {},
    })
    extract = AsyncMock(
        return_value=([{"gene_symbol": "TP53", "confidence": 0.9}], [{}])
    )
    merge = AsyncMock(
        return_value=({"gene_symbol": "TP53", "pathways": []}, {"input_tokens": 1})
    )
    # enrich returns the merged record unchanged.
    enrich = AsyncMock(side_effect=lambda gene, merged, cfg: merged)
    monkeypatch.setattr(pipeline, "run_fetch_stage", fetch)
    monkeypatch.setattr(pipeline, "run_extract_stage", extract)
    monkeypatch.setattr(pipeline, "run_merge_stage", merge)
    monkeypatch.setattr(pipeline, "run_enrich_stage", enrich)
    return fetch, extract, merge, enrich


@pytest.mark.asyncio
async def test_final_cache_hit_short_circuits(monkeypatch, config, stats):
    monkeypatch.setattr(
        pipeline, "read_cache",
        Mock(return_value={"gene_symbol": "TP53", "cached": True}),
    )
    read_raw = Mock(return_value=[{"x": 1}])
    monkeypatch.setattr(pipeline, "read_raw_cache", read_raw)
    fetch, extract, merge, enrich = _patch_stages(monkeypatch)

    result = await run_gene("TP53", config, stats)

    assert result["_from_cache"] is True
    assert stats.genes_cached == 1
    assert stats.genes_remerged == 0
    # A final-cache hit must not even consult the raw layer or any stage.
    read_raw.assert_not_called()
    fetch.assert_not_awaited()
    extract.assert_not_awaited()
    merge.assert_not_awaited()
    enrich.assert_not_awaited()


@pytest.mark.asyncio
async def test_raw_cache_hit_replays_merge_enrich(monkeypatch, config, stats, tmp_path):
    monkeypatch.setattr(pipeline, "ANNOTATIONS_JSONL", tmp_path / "ann.jsonl")
    monkeypatch.setattr(pipeline, "read_cache", Mock(return_value=None))
    monkeypatch.setattr(
        pipeline, "read_raw_cache",
        Mock(return_value=[{"gene_symbol": "TP53", "confidence": 0.9}]),
    )
    write_raw = Mock()
    write_final = Mock()
    monkeypatch.setattr(pipeline, "write_raw_cache", write_raw)
    monkeypatch.setattr(pipeline, "write_cache", write_final)
    fetch, extract, merge, enrich = _patch_stages(monkeypatch)

    result = await run_gene("TP53", config, stats)

    assert result == {"gene_symbol": "TP53", "pathways": []}
    assert stats.genes_remerged == 1
    assert stats.genes_cached == 0
    # Zero fetch/extract API work, and the raw cache is left untouched…
    fetch.assert_not_awaited()
    extract.assert_not_awaited()
    write_raw.assert_not_called()
    # …while merge + enrich are replayed and the final cache is rewritten.
    merge.assert_awaited_once()
    enrich.assert_awaited_once()
    write_final.assert_called_once()


@pytest.mark.asyncio
async def test_full_miss_runs_whole_chain(monkeypatch, config, stats, tmp_path):
    monkeypatch.setattr(pipeline, "ANNOTATIONS_JSONL", tmp_path / "ann.jsonl")
    monkeypatch.setattr(pipeline, "read_cache", Mock(return_value=None))
    monkeypatch.setattr(pipeline, "read_raw_cache", Mock(return_value=None))
    write_raw = Mock()
    write_final = Mock()
    monkeypatch.setattr(pipeline, "write_raw_cache", write_raw)
    monkeypatch.setattr(pipeline, "write_cache", write_final)
    fetch, extract, merge, enrich = _patch_stages(monkeypatch)

    result = await run_gene("TP53", config, stats)

    assert result == {"gene_symbol": "TP53", "pathways": []}
    assert stats.genes_cached == 0
    assert stats.genes_remerged == 0
    # The whole chain runs, the raw cache is written, then the final cache.
    fetch.assert_awaited_once()
    extract.assert_awaited_once()
    write_raw.assert_called_once()
    merge.assert_awaited_once()
    enrich.assert_awaited_once()
    write_final.assert_called_once()


@pytest.mark.asyncio
async def test_no_high_confidence_source_returns_none(monkeypatch, config, stats, tmp_path):
    # Merge returns (None, None) when no source cleared the gate → gene skipped,
    # nothing written to the final cache, but the (empty) raw cache is still saved.
    monkeypatch.setattr(pipeline, "ANNOTATIONS_JSONL", tmp_path / "ann.jsonl")
    monkeypatch.setattr(pipeline, "read_cache", Mock(return_value=None))
    monkeypatch.setattr(pipeline, "read_raw_cache", Mock(return_value=None))
    write_raw = Mock()
    write_final = Mock()
    monkeypatch.setattr(pipeline, "write_raw_cache", write_raw)
    monkeypatch.setattr(pipeline, "write_cache", write_final)
    _patch_stages(monkeypatch)
    monkeypatch.setattr(pipeline, "run_extract_stage", AsyncMock(return_value=([], [])))
    monkeypatch.setattr(pipeline, "run_merge_stage", AsyncMock(return_value=(None, None)))

    result = await run_gene("TP53", config, stats)

    assert result is None
    write_raw.assert_called_once()   # empty extractions are still cached
    write_final.assert_not_called()  # no final record to cache
