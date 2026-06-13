"""Centralized pipeline configuration.

A single dataclass that reads every environment variable used across the
pipeline, so configuration lives in one place rather than being scattered
through os.getenv() calls in each module. Import the shared ``config`` instance:

    from src.config import config

Scoring weights are validated to sum to 1.0 on construction.
"""

from dataclasses import dataclass, field
import os


@dataclass
class PipelineConfig:
    # Models
    extraction_model: str = os.getenv("EXTRACTION_MODEL", "claude-opus-4-8")
    merge_model: str = os.getenv("MERGE_MODEL", "claude-sonnet-4-6")
    max_tokens: int = int(os.getenv("MAX_TOKENS", "4096"))

    # Pipeline
    confidence_threshold: float = float(os.getenv("CONFIDENCE_THRESHOLD", "0.65"))
    pubmed_max_results: int = int(os.getenv("PUBMED_MAX_RESULTS", "20"))
    semaphore_limit: int = int(os.getenv("SEMAPHORE_LIMIT", "3"))
    log_level: str = os.getenv("LOG_LEVEL", "INFO")

    # Disease context
    disease_context: str = os.getenv("DISEASE_CONTEXT", "cancer")
    disease_terms: set = field(default_factory=lambda: set(
        os.getenv("DISEASE_TERMS",
        "cancer,tumor,carcinoma,sarcoma,lymphoma,leukemia,"
        "melanoma,adenocarcinoma,glioma,myeloma,blastoma,"
        "tumour,neoplasm,malignancy").split(",")
    ))

    # NCBI
    ncbi_email: str = os.getenv("NCBI_EMAIL", "")
    ncbi_api_key: str = os.getenv("NCBI_API_KEY", "")

    # STRING
    string_min_score: int = int(os.getenv("STRING_MIN_SCORE", "700"))
    string_limit: int = int(os.getenv("STRING_LIMIT", "50"))

    # GTEx
    gtex_tpm_threshold: float = float(os.getenv("GTEX_TPM_THRESHOLD", "10.0"))
    gtex_min_tissues: int = int(os.getenv("GTEX_MIN_TISSUES", "3"))

    # CellxGene
    enable_cellxgene: bool = os.getenv("ENABLE_CELLXGENE", "true").lower() == "true"
    census_version: str = os.getenv("CENSUS_VERSION", "2024-07-01")
    census_tissue: str = os.getenv("CENSUS_TISSUE", "lung")
    census_min_cells: int = int(os.getenv("CENSUS_MIN_CELLS", "50"))
    census_cache_dir: str = os.getenv("CENSUS_CACHE_DIR", "refs/census_cache/")

    # Scoring weights (must sum to 1.0)
    weight_betweenness: float = float(os.getenv("WEIGHT_BETWEENNESS", "0.25"))
    weight_degree: float = float(os.getenv("WEIGHT_DEGREE", "0.15"))
    weight_disease: float = float(os.getenv("WEIGHT_DISEASE", "0.35"))
    weight_druggability: float = float(os.getenv("WEIGHT_DRUGGABILITY", "0.10"))
    weight_cellxgene: float = float(os.getenv("WEIGHT_CELLXGENE", "0.15"))
    safety_penalty: float = float(os.getenv("SAFETY_PENALTY", "0.75"))

    # Fuzzy matching
    fuzzy_threshold: int = int(os.getenv("FUZZY_THRESHOLD", "85"))
    synonym_candidates: int = int(os.getenv("SYNONYM_CANDIDATES", "10"))
    synonym_model: str = os.getenv("SYNONYM_MODEL", "claude-sonnet-4-6")
    auto_update_synonyms: bool = os.getenv(
        "AUTO_UPDATE_SYNONYMS", "false").lower() == "true"

    # Layout
    layout_seed: int = int(os.getenv("LAYOUT_SEED", "42"))
    layout_k: float = float(os.getenv("LAYOUT_K", "2.5"))

    def __post_init__(self):
        total = (self.weight_betweenness + self.weight_degree +
                 self.weight_disease + self.weight_druggability +
                 self.weight_cellxgene)
        if abs(total - 1.0) > 0.001:
            raise ValueError(f"Scoring weights must sum to 1.0, got {total:.3f}")


config = PipelineConfig()
