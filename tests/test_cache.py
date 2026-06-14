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
import pipeline  # noqa: E402
from pipeline import (  # noqa: E402
    make_cache_key,
    make_extract_key,
    make_full_key,
    read_cache,
    write_cache,
    read_raw_cache,
    write_raw_cache,
    prune_stale_cache,
    _raw_cache_path,
    _final_cache_path,
)
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
    # safety_penalty is a scoring weight (applied in network.py), not an
    # extract/merge/enrich input, so it is part of NEITHER cache key.
    a = _config(safety_penalty=0.75)
    b = _config(safety_penalty=0.50)
    assert make_cache_key("TP53", a) == make_cache_key("TP53", b)
    assert make_extract_key("TP53", a) == make_extract_key("TP53", b)

def test_confidence_threshold_invalidates_final_not_raw():
    # The confidence gate is applied to the cached extractions at merge time, so a
    # threshold change must invalidate the final layer (re-filter + re-merge) but
    # NOT the raw layer (no re-extraction needed).
    a = _config(confidence_threshold=0.65)
    b = _config(confidence_threshold=0.90)
    assert make_full_key("TP53", a) != make_full_key("TP53", b)
    assert make_extract_key("TP53", a) == make_extract_key("TP53", b)


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


# --- Layer 1 (raw) key: scoped to fetch + extract inputs ONLY ---

def test_extract_key_excludes_merge_and_enrich_params():
    # The raw layer must NOT be invalidated by merge model, synonym/reference
    # files, or any enrich param — that is what makes a synonym edit free.
    ek = make_extract_key("TP53", _config())
    changed = _config(
        merge_model="claude-opus-4-8",
        census_tissue="liver",
        census_version="2099-01-01",
        string_min_score=400,
        gtex_tpm_threshold=99.0,
    )
    assert make_extract_key("TP53", changed) == ek

def test_extract_key_changes_with_extraction_inputs():
    ek = make_extract_key("TP53", _config())
    assert make_extract_key("BRCA1", _config()) != ek
    assert make_extract_key("TP53", _config(disease_context="fibrosis")) != ek
    assert make_extract_key("TP53", _config(extraction_model="x")) != ek
    assert make_extract_key("TP53", _config(pubmed_max_results=99)) != ek


# --- Layer 2 (full) key: sensitive to synonym/reference file content ---

def test_full_key_changes_with_synonym_file_content(tmp_path, monkeypatch):
    # Editing the synonym map must invalidate the final layer but NOT the raw layer.
    syn = tmp_path / "syn.json"
    ref = tmp_path / "ref.txt"
    syn.write_text("{}")
    ref.write_text("Signaling by WNT\n")
    monkeypatch.setattr(pipeline, "SYNONYMS_PATH", syn)
    monkeypatch.setattr(pipeline, "REACTOME_REF_PATH", str(ref))
    c = _config()

    fk_before = make_full_key("TP53", c)
    ek_before = make_extract_key("TP53", c)
    syn.write_text('{"foo bar": "Signaling by WNT"}')
    assert make_full_key("TP53", c) != fk_before     # final layer invalidated
    assert make_extract_key("TP53", c) == ek_before  # raw layer untouched

def test_full_key_changes_with_reactome_file_content(tmp_path, monkeypatch):
    syn = tmp_path / "syn.json"
    ref = tmp_path / "ref.txt"
    syn.write_text("{}")
    ref.write_text("Signaling by WNT\n")
    monkeypatch.setattr(pipeline, "SYNONYMS_PATH", syn)
    monkeypatch.setattr(pipeline, "REACTOME_REF_PATH", str(ref))
    c = _config()

    fk_before = make_full_key("TP53", c)
    ref.write_text("Signaling by WNT\nNew Canonical Pathway\n")
    assert make_full_key("TP53", c) != fk_before

def test_make_cache_key_is_full_key_alias():
    # Backward-compatible alias: the historical single key is the final-layer key.
    assert make_cache_key is make_full_key


# --- Layer 1 (raw) read/write: round-trip, subdirs, empty-list edge case ---

def test_raw_cache_round_trip_in_separate_subdir(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    extractions = [{"gene_symbol": "TP53", "confidence": 0.9}]
    write_raw_cache("TP53", c, extractions)
    write_cache("TP53", c, {"final": True})
    assert read_raw_cache("TP53", c) == extractions
    # The two layers live in distinct subdirectories.
    assert (tmp_path / "raw").is_dir()
    assert (tmp_path / "final").is_dir()

def test_raw_cache_empty_list_is_a_hit(tmp_path):
    # A gene with no high-confidence source caches []; that must read back as a hit
    # (return []), not as a miss (None) — so it isn't re-fetched on the next run.
    c = _config(cache_dir=str(tmp_path))
    write_raw_cache("TP53", c, [])
    assert read_raw_cache("TP53", c) == []

def test_enable_cache_false_disables_raw(tmp_path):
    c = _config(cache_dir=str(tmp_path), enable_cache=False)
    write_raw_cache("TP53", c, [{"a": 1}])  # no-op
    assert read_raw_cache("TP53", c) is None


# --- FORCE_REMERGE / FORCE_RERUN across the two layers ---

def test_force_remerge_misses_final_keeps_raw(tmp_path):
    # The free re-merge path: final cache bypassed, raw cache still served.
    c = _config(cache_dir=str(tmp_path))
    write_cache("TP53", c, {"final": True})
    write_raw_cache("TP53", c, [{"a": 1}])
    c.force_remerge = True
    assert read_cache("TP53", c) is None
    assert read_raw_cache("TP53", c) == [{"a": 1}]

def test_force_rerun_bypasses_both_layers(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    write_cache("TP53", c, {"final": True})
    write_raw_cache("TP53", c, [{"a": 1}])
    c.force_rerun = True
    assert read_cache("TP53", c) is None
    assert read_raw_cache("TP53", c) is None


# --- prune_stale_cache: delete old-key files, keep current/other-gene files ---

def _seed_current(gene, c):
    """Write the current-key raw + final files for a gene; return their names."""
    write_raw_cache(gene, c, [{"a": 1}])
    write_cache(gene, c, {"final": True})
    return _raw_cache_path(gene, c).name, _final_cache_path(gene, c).name

def test_prune_deletes_stale_keeps_current(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    raw, final = tmp_path / "raw", tmp_path / "final"
    cur_raw, cur_final = _seed_current("TP53", c)
    # Stale files (keys that don't match the current ones).
    (raw / "TP53_staaaaaaaaaa.json").write_text("{}")
    (final / "TP53_stbbbbbbbbbb.json").write_text("{}")

    raw_pruned, final_pruned = prune_stale_cache(c, ["TP53"])

    assert (raw_pruned, final_pruned) == (1, 1)
    assert (raw / cur_raw).exists()                       # current kept
    assert (final / cur_final).exists()
    assert not (raw / "TP53_staaaaaaaaaa.json").exists()  # stale deleted
    assert not (final / "TP53_stbbbbbbbbbb.json").exists()

def test_prune_leaves_other_genes_untouched(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    raw, final = tmp_path / "raw", tmp_path / "final"
    raw.mkdir(parents=True); final.mkdir(parents=True)
    # Files for a gene NOT in the prune list must survive (even with stale keys).
    (raw / "BRCA1_xxxxxxxxxxxx.json").write_text("{}")
    (final / "BRCA1_yyyyyyyyyyyy.json").write_text("{}")

    raw_pruned, final_pruned = prune_stale_cache(c, ["TP53"])

    assert (raw_pruned, final_pruned) == (0, 0)
    assert (raw / "BRCA1_xxxxxxxxxxxx.json").exists()
    assert (final / "BRCA1_yyyyyyyyyyyy.json").exists()

def test_prune_respects_prune_cache_false(tmp_path):
    c = _config(cache_dir=str(tmp_path), prune_cache=False)
    raw, final = tmp_path / "raw", tmp_path / "final"
    raw.mkdir(parents=True); final.mkdir(parents=True)
    (raw / "TP53_stale0000000.json").write_text("{}")
    (final / "TP53_stale1111111.json").write_text("{}")

    assert prune_stale_cache(c, ["TP53"]) == (0, 0)
    assert (raw / "TP53_stale0000000.json").exists()
    assert (final / "TP53_stale1111111.json").exists()

def test_prune_noop_when_cache_disabled(tmp_path):
    c = _config(cache_dir=str(tmp_path), enable_cache=False)
    raw = tmp_path / "raw"; raw.mkdir(parents=True)
    (raw / "TP53_stale0000000.json").write_text("{}")
    assert prune_stale_cache(c, ["TP53"]) == (0, 0)
    assert (raw / "TP53_stale0000000.json").exists()

def test_prune_returns_correct_counts(tmp_path):
    c = _config(cache_dir=str(tmp_path))
    raw, final = tmp_path / "raw", tmp_path / "final"
    _seed_current("TP53", c)  # current files (must be kept)
    for k in ("aaaaa", "bbbbb", "ccccc"):
        (raw / f"TP53_{k}.json").write_text("{}")    # 3 stale raw
    for k in ("ddddd", "eeeee"):
        (final / f"TP53_{k}.json").write_text("{}")  # 2 stale final

    assert prune_stale_cache(c, ["TP53"]) == (3, 2)
