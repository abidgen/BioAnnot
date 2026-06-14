"""Tests for RunStats: usage accounting, cost/cache math, and report rendering.

Importing ``pipeline`` runs ``load_dotenv()`` at module import, which would leak
.env values into ``os.environ`` for the rest of the session and break
``test_config.py``'s clean-env assertions. Snapshot ``os.environ``, import, then
restore so this module's import has no global side effect.
"""

import os as _os
from time import perf_counter

_ENV_SNAPSHOT = dict(_os.environ)
from pipeline import RunStats, _BOX_WIDTH  # noqa: E402
_os.environ.clear()
_os.environ.update(_ENV_SNAPSHOT)


# --- add_usage accumulation ---

def test_add_usage_accumulates():
    s = RunStats()
    s.add_usage({"input_tokens": 100, "output_tokens": 50,
                 "cache_read_input_tokens": 1000, "cache_creation_input_tokens": 200})
    s.add_usage({"input_tokens": 10, "output_tokens": 5,
                 "cache_read_input_tokens": 100, "cache_creation_input_tokens": 0})
    assert s.llm_calls == 2
    assert s.input_tokens == 110
    assert s.output_tokens == 55
    assert s.cache_read_tokens == 1100
    assert s.cache_created_tokens == 200


def test_add_usage_missing_keys_default_zero():
    s = RunStats()
    s.add_usage({})  # usage dict with no token keys
    assert s.llm_calls == 1
    assert s.input_tokens == 0
    assert s.output_tokens == 0
    assert s.cache_read_tokens == 0
    assert s.cache_created_tokens == 0


# --- cache_hit_rate ---

def test_cache_hit_rate_zero_calls():
    # No LLM calls → guard returns 0.0 (avoids div-by-zero).
    assert RunStats().cache_hit_rate == 0.0


def test_cache_hit_rate_computation():
    s = RunStats()
    s.llm_calls = 1
    s.cache_read_tokens = 900
    s.cache_created_tokens = 100
    assert s.cache_hit_rate == 0.9


# --- estimated_cost_usd (Sonnet pricing: $3/$15 in/out, $0.30/$3.75 cache r/w) ---

def test_estimated_cost_zero_when_empty():
    assert RunStats().estimated_cost_usd == 0.0


def test_estimated_cost_pricing():
    s = RunStats()
    s.input_tokens = 1_000_000
    s.output_tokens = 1_000_000
    s.cache_read_tokens = 1_000_000
    s.cache_created_tokens = 1_000_000
    # 3.0 + 15.0 + 0.30 + 3.75
    assert abs(s.estimated_cost_usd - 22.05) < 1e-9


# --- runtime_seconds ---

def test_runtime_seconds_nonnegative():
    s = RunStats()
    s.start_time = perf_counter() - 5.0
    rt = s.runtime_seconds
    assert isinstance(rt, float)
    assert rt >= 5.0


# --- print_report ---

def test_print_report_alignment(capsys):
    s = RunStats()
    s.genes_total = 5
    s.genes_succeeded = 5
    s.genes_cached = 2
    s.add_usage({"input_tokens": 45230, "output_tokens": 12840,
                 "cache_read_input_tokens": 22820, "cache_creation_input_tokens": 2697})
    anns = {"DEMO": {"pathways": [f"P{i}" for i in range(94)]
                      + [f"NON-CANONICAL: x{i}" for i in range(4)]}}
    s.print_report(anns)

    out = capsys.readouterr().out
    framed = [ln for ln in out.splitlines() if ln and ln[0] in "║╔╠╚"]
    assert framed, "no framed lines were rendered"
    widths = {len(ln) for ln in framed}
    assert widths == {_BOX_WIDTH + 2}, f"misaligned line widths: {widths}"


def test_print_report_counts_noncanonical(capsys):
    s = RunStats()
    anns = {
        "TP53": {"pathways": ["Signaling by WNT", "NON-CANONICAL: foo axis"]},
        "BRCA1": {"pathways": ["Signaling by BRCA1", "X", "NON-CANONICAL: bar"]},
    }
    s.print_report(anns)
    assert s.noncanonical_total == 2
    assert s.noncanonical_flagged == ["TP53: foo axis", "BRCA1: bar"]
