import pytest
import networkx as nx
from unittest.mock import patch
from src.config import PipelineConfig


def _safety(tier1=False, tier2=False):
    """Build a safety_assessment dict mirroring gtex_safety.assess_safety output."""
    if tier1:
        penalty = 0.60
    elif tier2:
        penalty = 0.80
    else:
        penalty = 1.00
    return {
        "safety_flag": tier1 or tier2,
        "tier1_flag": tier1,
        "tier2_flag": tier2,
        "tier1_high_tissues": {"Liver": 6.0} if tier1 else {},
        "tier2_high_tissues": {"Spleen": 15.0, "Whole Blood": 12.0} if tier2 else {},
        "high_expression_tissues": [],
        "max_vital_tpm": 6.0 if tier1 else 0.0,
        "max_tpm": 15.0 if (tier1 or tier2) else 0.0,
        "safety_penalty": penalty,
    }


def make_test_graph(safety_tp53=None, safety_brca1=None):
    safety_tp53 = safety_tp53 if safety_tp53 is not None else _safety()
    safety_brca1 = safety_brca1 if safety_brca1 is not None else _safety()
    G = nx.MultiDiGraph()
    G.add_node("TP53", node_type="target", confidence=1.0,
               disease_associations=[
                   {"disease": "cancer", "role": "tumor_suppressor",
                    "evidence_strength": "strong"}],
               druggability_notes="MDM2 inhibitors",
               cellxgene_expression={"cell_type_count": 5},
               safety_assessment=safety_tp53)
    G.add_node("BRCA1", node_type="target", confidence=1.0,
               disease_associations=[
                   {"disease": "cancer", "role": "tumor_suppressor",
                    "evidence_strength": "strong"}],
               druggability_notes="PARP inhibitors",
               cellxgene_expression={"cell_type_count": 5},
               safety_assessment=safety_brca1)
    return G


_WEIGHTS = {
    "WEIGHT_BETWEENNESS": "0.25", "WEIGHT_DEGREE": "0.15",
    "WEIGHT_DISEASE": "0.35", "WEIGHT_DRUGGABILITY": "0.10",
    "WEIGHT_CELLXGENE": "0.15",
}


def _set_weights(monkeypatch):
    for k, v in _WEIGHTS.items():
        monkeypatch.setenv(k, v)


def test_safety_penalty_reduces_score(monkeypatch):
    _set_weights(monkeypatch)
    from src.network import compute_priority_scores
    # TP53 flagged (tier1), BRCA1 clean — identical otherwise.
    G = make_test_graph(safety_tp53=_safety(tier1=True), safety_brca1=_safety())
    scores = compute_priority_scores(G)
    tp53 = next(s for s in scores if s["gene"] == "TP53")
    brca1 = next(s for s in scores if s["gene"] == "BRCA1")
    assert tp53["composite"] < brca1["composite"]
    assert tp53["safety_flag"] is True
    assert brca1["safety_flag"] is False


def test_no_safety_penalty_when_not_flagged(monkeypatch):
    _set_weights(monkeypatch)
    from src.network import compute_priority_scores
    G = make_test_graph(safety_tp53=_safety(), safety_brca1=_safety())
    scores = compute_priority_scores(G)
    for s in scores:
        assert s["safety_penalty_applied"] is False
        assert s["safety_penalty"] == 1.0


def test_tier1_penalty_is_060(monkeypatch):
    _set_weights(monkeypatch)
    from src.network import compute_priority_scores
    G = make_test_graph(safety_tp53=_safety(tier1=True), safety_brca1=_safety())
    scores = compute_priority_scores(G)
    tp53 = next(s for s in scores if s["gene"] == "TP53")
    brca1 = next(s for s in scores if s["gene"] == "BRCA1")
    assert tp53["tier1_flag"] is True
    assert tp53["safety_penalty"] == 0.60
    # Identical inputs but for the penalty → exact 0.60 ratio.
    assert tp53["composite"] == pytest.approx(brca1["composite"] * 0.60)


def test_tier2_only_penalty_is_080(monkeypatch):
    _set_weights(monkeypatch)
    from src.network import compute_priority_scores
    G = make_test_graph(safety_tp53=_safety(tier2=True), safety_brca1=_safety())
    scores = compute_priority_scores(G)
    tp53 = next(s for s in scores if s["gene"] == "TP53")
    brca1 = next(s for s in scores if s["gene"] == "BRCA1")
    assert tp53["tier1_flag"] is False
    assert tp53["tier2_flag"] is True
    assert tp53["safety_penalty"] == 0.80
    assert tp53["composite"] == pytest.approx(brca1["composite"] * 0.80)


def test_tier1_precedence_over_tier2(monkeypatch):
    _set_weights(monkeypatch)
    from src.network import compute_priority_scores
    # Both tiers fire → tier 1 (0.60) wins.
    G = make_test_graph(safety_tp53=_safety(tier1=True, tier2=True))
    scores = compute_priority_scores(G)
    tp53 = next(s for s in scores if s["gene"] == "TP53")
    assert tp53["tier1_flag"] is True
    assert tp53["tier2_flag"] is True
    assert tp53["safety_penalty"] == 0.60


def test_only_target_nodes_scored():
    from src.network import compute_priority_scores
    G = make_test_graph()
    G.add_node("MDM2", node_type="interactor", confidence=None)
    scores = compute_priority_scores(G)
    genes = [s["gene"] for s in scores]
    assert "MDM2" not in genes
    assert "TP53" in genes
    assert "BRCA1" in genes
