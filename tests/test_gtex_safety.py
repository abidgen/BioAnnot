"""Unit tests for the two-tier GTEx normal-tissue safety filter.

These inject a synthetic median-TPM frame into the module-level ``_GTEX`` cache so
the tier logic is tested without downloading the real GTEx table.
"""

import pandas as pd
import pytest

import src.filters.gtex_safety as gs


def _fake_gtex(rows: dict[str, dict[str, float]]) -> pd.DataFrame:
    """Build a Description-indexed frame; unspecified sensitive tissues → 0.0 TPM."""
    tissues = sorted(gs.SENSITIVE_TISSUES)
    data = {
        gene: {t: tpms.get(t, 0.0) for t in tissues}
        for gene, tpms in rows.items()
    }
    frame = pd.DataFrame.from_dict(data, orient="index")
    frame.index.name = "Description"
    return frame


@pytest.fixture
def gtex(monkeypatch):
    """Return a helper that installs a fake _GTEX table for the test."""
    def _install(rows):
        monkeypatch.setattr(gs, "_GTEX", _fake_gtex(rows))
    return _install


def test_tier1_flag_at_any_vital_organ_over_5(gtex):
    # Single vital organ (Liver) just over 5 TPM → hard flag.
    gtex({"GENE1": {"Liver": 6.0}})
    r = gs.assess_safety("GENE1")
    assert r["tier1_flag"] is True
    assert r["tier2_flag"] is False
    assert r["safety_flag"] is True
    assert r["safety_penalty"] == 0.60
    assert r["tier1_high_tissues"] == {"Liver": 6.0}
    assert r["max_vital_tpm"] == 6.0


def test_tier1_threshold_is_strict(gtex):
    # Exactly 5.0 is NOT > 5.0 → no flag.
    gtex({"EDGE": {"Liver": 5.0}})
    r = gs.assess_safety("EDGE")
    assert r["tier1_flag"] is False
    assert r["safety_penalty"] == 1.00


def test_tier2_flag_at_two_secondary_tissues_over_10(gtex):
    gtex({"GENE2": {"Spleen": 15.0, "Whole Blood": 12.0}})
    r = gs.assess_safety("GENE2")
    assert r["tier1_flag"] is False
    assert r["tier2_flag"] is True
    assert r["safety_flag"] is True
    assert r["safety_penalty"] == 0.80
    assert set(r["tier2_high_tissues"]) == {"Spleen", "Whole Blood"}


def test_tier2_needs_at_least_two_tissues(gtex):
    # Only one secondary tissue over 10 → below tier2_min_tissues (2) → no flag.
    gtex({"GENE3": {"Spleen": 15.0}})
    r = gs.assess_safety("GENE3")
    assert r["tier2_flag"] is False
    assert r["safety_flag"] is False
    assert r["safety_penalty"] == 1.00


def test_tier1_takes_precedence_over_tier2(gtex):
    # Both tiers fire simultaneously → tier 1 penalty (0.60) wins.
    gtex({"GENE4": {"Liver": 6.0, "Spleen": 15.0, "Whole Blood": 12.0}})
    r = gs.assess_safety("GENE4")
    assert r["tier1_flag"] is True
    assert r["tier2_flag"] is True
    assert r["safety_penalty"] == 0.60


def test_gene_absent_from_gtex_is_unflagged(gtex):
    gtex({"SOMETHING": {"Liver": 100.0}})
    r = gs.assess_safety("NOT_IN_TABLE")
    assert r["safety_flag"] is False
    assert r["tier1_flag"] is False
    assert r["tier2_flag"] is False
    assert r["safety_penalty"] == 1.00
    assert r["max_tpm"] == 0.0
    assert r["max_vital_tpm"] == 0.0


def test_table_unavailable_degrades_to_no_concern(monkeypatch):
    monkeypatch.setattr(gs, "_GTEX", None)
    r = gs.assess_safety("ANY")
    assert r["safety_flag"] is False
    assert r["safety_penalty"] == 1.00


def test_penalty_values_are_config_driven(gtex, monkeypatch):
    # Env-override both penalties; a fresh config picks them up and they flow
    # through assess_safety (tier 1 → tier1_penalty, tier 2-only → tier2_penalty).
    monkeypatch.setenv("GTEX_TIER1_PENALTY", "0.4")
    monkeypatch.setenv("GTEX_TIER2_PENALTY", "0.7")
    from src.config import PipelineConfig
    cfg = PipelineConfig()
    assert cfg.gtex_tier1_penalty == 0.4
    assert cfg.gtex_tier2_penalty == 0.7

    gtex({
        "VITAL": {"Liver": 6.0},
        "SECONDARY": {"Spleen": 15.0, "Whole Blood": 12.0},
    })
    r1 = gs.assess_safety(
        "VITAL",
        tier1_penalty=cfg.gtex_tier1_penalty,
        tier2_penalty=cfg.gtex_tier2_penalty,
    )
    r2 = gs.assess_safety(
        "SECONDARY",
        tier1_penalty=cfg.gtex_tier1_penalty,
        tier2_penalty=cfg.gtex_tier2_penalty,
    )
    assert r1["tier1_flag"] is True and r1["safety_penalty"] == 0.4
    assert r2["tier2_flag"] is True and r2["safety_penalty"] == 0.7


def test_high_expression_tissues_is_union_of_both_tiers(gtex):
    gtex({"GENE5": {"Liver": 6.0, "Spleen": 15.0, "Whole Blood": 12.0}})
    r = gs.assess_safety("GENE5")
    assert r["high_expression_tissues"] == sorted(
        ["Liver", "Spleen", "Whole Blood"]
    )
    assert r["max_tpm"] == 15.0
