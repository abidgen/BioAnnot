"""Tests for run_gene_with_retry: retry on transient errors, fail-fast otherwise.

run_gene is replaced with an AsyncMock and asyncio.sleep is patched so the 2s/4s
backoff doesn't actually wait. Importing ``pipeline`` runs ``load_dotenv()`` at
import; snapshot/restore ``os.environ`` so it has no global side effect (keeps
``test_config.py`` clean-env assertions valid).
"""

import os as _os
from unittest.mock import AsyncMock

import httpx
import pytest

_ENV_SNAPSHOT = dict(_os.environ)
import pipeline  # noqa: E402
from pipeline import run_gene_with_retry, RunStats  # noqa: E402
_os.environ.clear()
_os.environ.update(_ENV_SNAPSHOT)

from src.config import PipelineConfig  # noqa: E402


@pytest.fixture
def config():
    return PipelineConfig()


@pytest.fixture
def stats():
    return RunStats()


@pytest.mark.asyncio
async def test_retry_succeeds_first_attempt(monkeypatch, config, stats):
    mock = AsyncMock(return_value={"gene_symbol": "TP53"})
    monkeypatch.setattr(pipeline, "run_gene", mock)

    result = await run_gene_with_retry("TP53", config, stats)

    assert result == {"gene_symbol": "TP53"}
    assert mock.await_count == 1


@pytest.mark.asyncio
async def test_retry_recovers_after_transient(monkeypatch, config, stats):
    mock = AsyncMock(
        side_effect=[httpx.ConnectTimeout("boom"), {"gene_symbol": "TP53"}]
    )
    monkeypatch.setattr(pipeline, "run_gene", mock)
    sleep_mock = AsyncMock()
    monkeypatch.setattr("asyncio.sleep", sleep_mock)

    result = await run_gene_with_retry("TP53", config, stats)

    assert result == {"gene_symbol": "TP53"}
    assert mock.await_count == 2
    sleep_mock.assert_awaited_once_with(2)  # 2 ** 1


@pytest.mark.asyncio
async def test_retry_exhausts_and_raises(monkeypatch, config, stats):
    mock = AsyncMock(side_effect=httpx.ConnectError("down"))
    monkeypatch.setattr(pipeline, "run_gene", mock)
    sleep_mock = AsyncMock()
    monkeypatch.setattr("asyncio.sleep", sleep_mock)

    with pytest.raises(httpx.ConnectError):
        await run_gene_with_retry("TP53", config, stats)

    assert mock.await_count == 3  # all 3 attempts used
    # Backoff after attempts 1 and 2, exponentially: 2s then 4s.
    assert [c.args[0] for c in sleep_mock.await_args_list] == [2, 4]


@pytest.mark.asyncio
async def test_non_retryable_fails_fast(monkeypatch, config, stats):
    mock = AsyncMock(side_effect=ValueError("bad data"))
    monkeypatch.setattr(pipeline, "run_gene", mock)
    sleep_mock = AsyncMock()
    monkeypatch.setattr("asyncio.sleep", sleep_mock)

    with pytest.raises(ValueError):
        await run_gene_with_retry("TP53", config, stats)

    assert mock.await_count == 1  # no retry on a non-retryable error
    sleep_mock.assert_not_awaited()
