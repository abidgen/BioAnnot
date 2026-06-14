"""Centralized pipeline configuration.

A single dataclass that reads every environment variable used across the
pipeline, so configuration lives in one place rather than being scattered
through os.getenv() calls in each module. Import the shared ``config`` instance:

    from src.config import config

Each field reads its environment variable through ``default_factory``, so
``PipelineConfig()`` re-reads the current environment on every instantiation
(rather than freezing values at class-definition time). Scoring weights are
validated to sum to 1.0 on construction.
"""

from dataclasses import dataclass, field
import os


@dataclass
class PipelineConfig:
    # Models
    extraction_model: str = field(
        default_factory=lambda: os.getenv("EXTRACTION_MODEL", "claude-opus-4-8"))
    merge_model: str = field(
        default_factory=lambda: os.getenv("MERGE_MODEL", "claude-sonnet-4-6"))
    max_tokens: int = field(
        default_factory=lambda: int(os.getenv("MAX_TOKENS", "4096")))

    # Pipeline
    confidence_threshold: float = field(
        default_factory=lambda: float(os.getenv("CONFIDENCE_THRESHOLD", "0.65")))
    pubmed_max_results: int = field(
        default_factory=lambda: int(os.getenv("PUBMED_MAX_RESULTS", "20")))
    semaphore_limit: int = field(
        default_factory=lambda: int(os.getenv("SEMAPHORE_LIMIT", "3")))
    log_level: str = field(
        default_factory=lambda: os.getenv("LOG_LEVEL", "INFO"))

    # Disease context
    disease_context: str = field(
        default_factory=lambda: os.getenv("DISEASE_CONTEXT", "cancer"))
    disease_terms: set = field(default_factory=lambda: set(
        os.getenv("DISEASE_TERMS",
        "cancer,tumor,carcinoma,sarcoma,lymphoma,leukemia,"
        "melanoma,adenocarcinoma,glioma,myeloma,blastoma,"
        "tumour,neoplasm,malignancy").split(",")
    ))

    # NCBI
    ncbi_email: str = field(
        default_factory=lambda: os.getenv("NCBI_EMAIL", ""))
    ncbi_api_key: str = field(
        default_factory=lambda: os.getenv("NCBI_API_KEY", ""))

    # STRING
    string_min_score: int = field(
        default_factory=lambda: int(os.getenv("STRING_MIN_SCORE", "700")))
    string_limit: int = field(
        default_factory=lambda: int(os.getenv("STRING_LIMIT", "50")))

    # GTEx
    gtex_tpm_threshold: float = field(
        default_factory=lambda: float(os.getenv("GTEX_TPM_THRESHOLD", "10.0")))
    gtex_min_tissues: int = field(
        default_factory=lambda: int(os.getenv("GTEX_MIN_TISSUES", "3")))

    # CellxGene
    enable_cellxgene: bool = field(
        default_factory=lambda: os.getenv("ENABLE_CELLXGENE", "true").lower() == "true")
    census_version: str = field(
        default_factory=lambda: os.getenv("CENSUS_VERSION", "2024-07-01"))
    census_tissue: str = field(
        default_factory=lambda: os.getenv("CENSUS_TISSUE", "lung"))
    census_min_cells: int = field(
        default_factory=lambda: int(os.getenv("CENSUS_MIN_CELLS", "50")))
    census_cache_dir: str = field(
        default_factory=lambda: os.getenv("CENSUS_CACHE_DIR", "refs/census_cache/"))

    # Scoring weights (must sum to 1.0)
    weight_betweenness: float = field(
        default_factory=lambda: float(os.getenv("WEIGHT_BETWEENNESS", "0.25")))
    weight_degree: float = field(
        default_factory=lambda: float(os.getenv("WEIGHT_DEGREE", "0.15")))
    weight_disease: float = field(
        default_factory=lambda: float(os.getenv("WEIGHT_DISEASE", "0.35")))
    weight_druggability: float = field(
        default_factory=lambda: float(os.getenv("WEIGHT_DRUGGABILITY", "0.10")))
    weight_cellxgene: float = field(
        default_factory=lambda: float(os.getenv("WEIGHT_CELLXGENE", "0.15")))
    safety_penalty: float = field(
        default_factory=lambda: float(os.getenv("SAFETY_PENALTY", "0.75")))

    # Fuzzy matching
    fuzzy_threshold: int = field(
        default_factory=lambda: int(os.getenv("FUZZY_THRESHOLD", "85")))
    synonym_candidates: int = field(
        default_factory=lambda: int(os.getenv("SYNONYM_CANDIDATES", "10")))
    synonym_model: str = field(
        default_factory=lambda: os.getenv("SYNONYM_MODEL", "claude-sonnet-4-6"))
    auto_update_synonyms: bool = field(
        default_factory=lambda: os.getenv("AUTO_UPDATE_SYNONYMS", "false").lower() == "true")

    # Layout
    layout_seed: int = field(
        default_factory=lambda: int(os.getenv("LAYOUT_SEED", "42")))
    layout_k: float = field(
        default_factory=lambda: float(os.getenv("LAYOUT_K", "2.5")))

    # On-disk resume cache (skip a gene's stages if a prior run cached it)
    enable_cache: bool = field(
        default_factory=lambda:
        os.getenv("ENABLE_CACHE", "true").lower() == "true"
    )
    cache_dir: str = field(
        default_factory=lambda: os.getenv("CACHE_DIR", "outputs/cache/")
    )
    # Force a fresh run: bypass BOTH cache layers (raw + final) and recompute the
    # whole chain, rewriting both caches.
    force_rerun: bool = field(
        default_factory=lambda:
        os.getenv("FORCE_RERUN", "false").lower() == "true"
    )
    # Force a re-merge: bypass the final cache only, replay merge + enrich from the
    # raw extraction cache, and rewrite the final cache. Zero API cost when the raw
    # cache is warm (no fetch/extract). Use after editing the synonym map or the
    # Reactome reference. Ignored when force_rerun is set (that bypasses everything).
    force_remerge: bool = field(
        default_factory=lambda:
        os.getenv("FORCE_REMERGE", "false").lower() == "true"
    )

    def __post_init__(self):
        total = (self.weight_betweenness + self.weight_degree +
                 self.weight_disease + self.weight_druggability +
                 self.weight_cellxgene)
        if abs(total - 1.0) > 0.001:
            raise ValueError(f"Scoring weights must sum to 1.0, got {total:.3f}")


config = PipelineConfig()
