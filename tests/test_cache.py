"""Cache invalidation tests for the on-disk resume cache (pipeline.py).

``make_cache_key`` digests the inputs that determine a gene's output, so changing
any of them must yield a different key (invalidating the cache), while changing a
field outside the key must not.

Importing ``pipeline`` runs ``load_dotenv()`` at module import, which would leak
.env values into ``os.environ`` for the rest of the test session and break
``test_config.py``'s clean-env assertions. Snapshot ``os.environ``, import, then
restore so this module's import has no global side effect.
"""

import os as _os

_ENV_SNAPSHOT = dict(_os.environ)
from pipeline import make_cache_key, read_cache, write_cache  # noqa: E402
_os.environ.clear()
_os.environ.update(_ENV_SNAPSHOT)

from src.config import PipelineConfig  # noqa: E402


def _config(**overrides):
    """A PipelineConfig with explicit field overrides (independent of ambient env)."""
    c = PipelineConfig()
    for key, value in overrides.items():
        setattr(c, key, value)
    return c


# --- make_cache_key: stability and per-component invalidation ---

def test_cache_key_stable_for_same_config():
    c = _config()
    assert make_cache_key("TP53", c) == make_cache_key("TP53", c)

def test_cache_key_differs_per_gene():
    c = _config()
    assert make_cache_key("TP53", c) != make_cache_key("BRCA1", c)

def test_cache_key_changes_with_extraction_model():
    a = _config(extraction_model="claude-opus-4-8")
    b = _config(extraction_model="claude-sonnet-4-6")
    assert make_cache_key("TP53", a) != make_cache_key("TP53", b)

def test_cache_key_changes_with_merge_model():
    a = _config(merge_model="claude-sonnet-4-6")
    b = _config(merge_model="claude-opus-4-8")
    assert make_cache_key("TP53", a) != make_cache_key("TP53", b)

def test_cache_key_changes_with_disease_context():
    a = _config(disease_context="cancer")
    b = _config(disease_context="fibrosis")
    assert make_cache_key("TP53", a) != make_cache_key("TP53", b)

def test_cache_key_changes_with_pubmed_max_results():
    a = _config(pubmed_max_results=20)
    b = _config(pubmed_max_results=50)
    assert make_cache_key("TP53", a) != make_cache_key("TP53", b)

def test_cache_key_changes_with_census_tissue():
    a = _config(census_tissue="lung")
    b = _config(census_tissue="liver")
    assert make_cache_key("TP53", a) != make_cache_key("TP53", b)

def test_cache_key_changes_with_census_version():
    a = _config(census_version="2024-07-01")
    b = _config(census_version="2023-12-15")
    assert make_cache_key("TP53", a) != make_cache_key("TP53", b)

def test_cache_key_unaffected_by_non_key_fields():
    # confidence_threshold and safety_penalty are NOT part of the cache key.
    a = _config(confidence_threshold=0.65, safety_penalty=0.75)
    b = _config(confidence_threshold=0.90, safety_penalty=0.50)
    assert make_cache_key("TP53", a) == make_cache_key("TP53", b)


# --- read_cache / write_cache: round-trip and invalidation on disk ---

def test_read_cache_hit_after_write(tmp_path):
    c = _config(cache_dir=str(tmp_path), extraction_model="claude-opus-4-8")
    result = {"gene_symbol": "TP53", "confidence": 0.9}
    write_cache("TP53", c, result)
    assert read_cache("TP53", c) == result

def test_read_cache_invalidated_when_config_changes(tmp_path):
    # Write under config A, then read under a config whose key-component changed:
    # the new key points at a different (nonexistent) file → cache miss.
    a = _config(cache_dir=str(tmp_path), extraction_model="claude-opus-4-8")
    write_cache("TP53", a, {"gene_symbol": "TP53", "confidence": 0.9})
    b = _config(cache_dir=str(tmp_path), extraction_model="claude-sonnet-4-6")
    assert read_cache("TP53", b) is None
    # The original (unchanged) config still hits.
    assert read_cache("TP53", a) is not None

def test_force_rerun_bypasses_cache(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    write_cache("TP53", c, {"x": 1})
    assert read_cache("TP53", c) is not None
    c.force_rerun = True
    assert read_cache("TP53", c) is None

def test_enable_cache_false_disables_read(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    write_cache("TP53", c, {"x": 1})
    c.enable_cache = False
    assert read_cache("TP53", c) is None
