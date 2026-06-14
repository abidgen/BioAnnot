import pytest
import networkx as nx
from unittest.mock import patch
from src.config import PipelineConfig

def make_test_graph(safety_flag_tp53=True, safety_flag_brca1=False):
    G = nx.MultiDiGraph()
    G.add_node("TP53", node_type="target", confidence=1.0,
               disease_associations=[
                   {"disease": "cancer", "role": "tumor_suppressor",
                    "evidence_strength": "strong"}],
               druggability_notes="MDM2 inhibitors",
               cellxgene_expression={"cell_type_count": 5},
               safety_assessment={"safety_flag": safety_flag_tp53})
    G.add_node("BRCA1", node_type="target", confidence=1.0,
               disease_associations=[
                   {"disease": "cancer", "role": "tumor_suppressor",
                    "evidence_strength": "strong"}],
               druggability_notes="PARP inhibitors",
               cellxgene_expression={"cell_type_count": 5},
               safety_assessment={"safety_flag": safety_flag_brca1})
    return G

def test_safety_penalty_reduces_score(monkeypatch):
    monkeypatch.setenv("WEIGHT_BETWEENNESS", "0.25")
    monkeypatch.setenv("WEIGHT_DEGREE", "0.15")
    monkeypatch.setenv("WEIGHT_DISEASE", "0.35")
    monkeypatch.setenv("WEIGHT_DRUGGABILITY", "0.10")
    monkeypatch.setenv("WEIGHT_CELLXGENE", "0.15")
    monkeypatch.setenv("SAFETY_PENALTY", "0.75")
    from src.network import compute_priority_scores
    G = make_test_graph(safety_flag_tp53=True, safety_flag_brca1=False)
    scores = compute_priority_scores(G)
    tp53 = next(s for s in scores if s["gene"] == "TP53")
    brca1 = next(s for s in scores if s["gene"] == "BRCA1")
    assert tp53["composite"] < brca1["composite"]
    assert tp53["safety_flag"] == True
    assert brca1["safety_flag"] == False

def test_no_safety_penalty_when_not_flagged(monkeypatch):
    monkeypatch.setenv("WEIGHT_BETWEENNESS", "0.25")
    monkeypatch.setenv("WEIGHT_DEGREE", "0.15")
    monkeypatch.setenv("WEIGHT_DISEASE", "0.35")
    monkeypatch.setenv("WEIGHT_DRUGGABILITY", "0.10")
    monkeypatch.setenv("WEIGHT_CELLXGENE", "0.15")
    monkeypatch.setenv("SAFETY_PENALTY", "0.75")
    from src.network import compute_priority_scores
    G = make_test_graph(safety_flag_tp53=False, safety_flag_brca1=False)
    scores = compute_priority_scores(G)
    for s in scores:
        assert s["safety_penalty_applied"] == False

def test_only_target_nodes_scored():
    from src.network import compute_priority_scores
    G = make_test_graph()
    G.add_node("MDM2", node_type="interactor", confidence=None)
    scores = compute_priority_scores(G)
    genes = [s["gene"] for s in scores]
    assert "MDM2" not in genes
    assert "TP53" in genes
    assert "BRCA1" in genes
