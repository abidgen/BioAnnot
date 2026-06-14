import os
import pytest
from src.config import PipelineConfig

def test_config_defaults():
    config = PipelineConfig()
    assert config.extraction_model == "claude-opus-4-8"
    assert config.merge_model == "claude-sonnet-4-6"
    assert config.fuzzy_threshold == 85
    assert config.safety_penalty == 0.75
    assert config.confidence_threshold == 0.65
    assert config.pubmed_max_results == 20
    assert config.enable_cellxgene == True

def test_config_env_override(monkeypatch):
    monkeypatch.setenv("FUZZY_THRESHOLD", "90")
    monkeypatch.setenv("EXTRACTION_MODEL", "claude-sonnet-4-6")
    monkeypatch.setenv("SAFETY_PENALTY", "0.5")
    config = PipelineConfig()
    assert config.fuzzy_threshold == 90
    assert config.extraction_model == "claude-sonnet-4-6"
    assert config.safety_penalty == 0.5

def test_config_weight_validation_passes():
    config = PipelineConfig()
    total = (config.weight_betweenness + config.weight_degree +
             config.weight_disease + config.weight_druggability +
             config.weight_cellxgene)
    assert abs(total - 1.0) < 0.001

def test_config_weight_validation_fails(monkeypatch):
    monkeypatch.setenv("WEIGHT_BETWEENNESS", "0.5")
    monkeypatch.setenv("WEIGHT_DEGREE", "0.5")
    monkeypatch.setenv("WEIGHT_DISEASE", "0.5")
    monkeypatch.setenv("WEIGHT_DRUGGABILITY", "0.0")
    monkeypatch.setenv("WEIGHT_CELLXGENE", "0.0")
    with pytest.raises(ValueError, match="must sum to 1.0"):
        PipelineConfig()

def test_config_disease_terms_default():
    config = PipelineConfig()
    assert "cancer" in config.disease_terms
    assert "leukemia" in config.disease_terms
    assert "carcinoma" in config.disease_terms

def test_config_enable_cellxgene_false(monkeypatch):
    monkeypatch.setenv("ENABLE_CELLXGENE", "false")
    config = PipelineConfig()
    assert config.enable_cellxgene == False
