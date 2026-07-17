# Pipeline Analysis & Internals

A deep reference for how the annotation pipeline actually works — the per-stage logic, the
scoring model, the caching machinery, and (importantly) the **caveats** that affect how the
outputs should be read. For setup and a quick tour see [`README.md`](../README.md); this
document is the "what it really does and where it bends" companion.

---

## 1. What the pipeline computes

Given a list of gene symbols, the pipeline produces, per gene, a single reconciled annotation
record (functions, cellular states, pathways, disease associations, interactors, druggability
notes, confidence) enriched with database-backed signals (STRING PPIs, a GTEx normal-tissue
safety flag, CellxGene single-cell expression), then builds a NetworkX graph and ranks the
genes by a composite prioritization score.

Two ideas run through everything:

- **LLM where judgement is needed, lookups where facts exist.** Extraction and merge use the
  Anthropic API with forced tool use; STRING / GTEx / CellxGene are deterministic lookups that
  bypass the LLM and attach directly to the record.
- **Disk is the source of truth.** Every stage writes to the run's timestamped directory,
  `outputs/runs/{timestamp}/` (with `outputs/latest` pointing at it). Re-runs resume from a
  two-layer cache — shared across runs under `outputs/cache/`, outside the per-run directory —
  rather than recomputing.

---

## 2. End-to-end flow

```
inputs/target_genes.txt   ─┐
refs/reactome_pathways.txt │  loaded once
refs/pathway_synonyms.json ┘
        │
        ▼  per gene, concurrency = SEMAPHORE_LIMIT (default 3), with retry + error isolation
┌────────────────────────────────────────────────────────────────────────┐
│ run_gene(gene)                                                         │
│   1. final cache (full_key)  ── HIT ▶ return record  (genes_cached)    │
│   2. raw cache  (extract_key)── HIT ▶ skip 3, go to 4 (genes_remerged) │
│   3. FETCH + EXTRACT ▶ write raw cache              (full chain)       │
│        PubMed → UniProt → OpenTargets ▶ LLM extract per source         │
│   4. confidence gate (applied here, not at extract time)               │
│   5. MERGE  (LLM for ≥2 sources; local for 1) ▶ pathway canonicalize   │
│   6. ENRICH  STRING + GTEx safety + CellxGene                          │
│   7. write final cache; append to annotations.jsonl  (genes_succeeded) │
└────────────────────────────────────────────────────────────────────────┘
        │  all genes complete
        ▼
  prune_stale_cache → final_annotations.json
        ▼
  build_target_network → compute_priority_scores → prioritized_targets.tsv + .gpickle
        ▼
  (optional) AUTO_UPDATE_SYNONYMS → run report
```

`pipeline.py` orchestrates this; `batch_pipeline.py` is a Batch-API variant for ≥50 genes
(see [§9](#9-known-gaps--differences)).

---

## 3. Stage logic

### 3.1 Fetch (`run_fetch_stage`)
- **PubMed**: ESearch for `"{gene}[gene] AND ({disease OR-clause})"` (built from
  `DISEASE_CONTEXT` + `DISEASE_TERMS`), up to `PUBMED_MAX_RESULTS` PMIDs, validated to 7–8
  digits, then EFetch abstracts. Missing abstracts are logged and skipped.
- **UniProt**: reviewed human entry (`gene_exact` + `organism_id:9606 + reviewed:true`) →
  function text, subcellular locations, GO terms, keywords.
- **OpenTargets**: GraphQL — resolve symbol → Ensembl ID, then function descriptions,
  pathways, top disease associations.

All fetchers are async with retry; a source with no entry returns `{}` and is simply dropped.

### 3.2 Extract (`run_extract_stage`)
One LLM call **per source** (`EXTRACTION_MODEL`, default `claude-opus-4-8`) at
**`temperature=0`** with the `annotate_target` tool **forced** via `tool_choice`, so output is
always structured. Input text is truncated to `EXTRACTION_MAX_WORDS` (default 5000). The system prompt embeds the disease
context. Every successful per-source extraction is returned **unfiltered** and persisted to
`{run}/raw/{gene}_raw.json` (in the run's timestamped directory) for provenance.

> PMIDs are only ever passed through from the validated fetcher list — the model is never
> trusted to produce citations.

### 3.3 Confidence gate (in `run_gene`, at consumption)
Sources below `CONFIDENCE_THRESHOLD` (default 0.65) are dropped **here**, not at extract time.
This is deliberate: the raw cache stores the *unfiltered* extractions, so changing the
threshold takes effect on a re-merge (it is part of `full_key`) without re-extracting. See
[§5](#5-two-layer-cache).

### 3.4 Merge (`merge_annotations`)
- **Single source** → returned directly (no LLM call), with pathway canonicalization applied.
- **Multiple sources** → one LLM call (`MERGE_MODEL`, default `claude-sonnet-4-6`) at
  **`temperature=0`** with the same forced tool, reconciling conflicts by evidence strength and
  source priority
  (UniProt > OpenTargets > PubMed for function; Reactome > OpenTargets > PubMed for pathways),
  unioning interactors, and noting context-dependence for conflicting disease roles.
- After the model returns, pathway names are canonicalized and `source_pmids` are **unioned
  from the input sources** (the tool schema has no PMID field, preventing fabricated cites).

### 3.5 Pathway canonicalization (`src/pathways.py`)
Each pathway name is resolved against the Reactome reference in priority order:
1. **exact** — normalized (lowercase, trimmed, Reactome `(R-HSA-…)` suffix stripped);
2. **synonym** — `refs/pathway_synonyms.json` (informal → canonical, validated on load);
3. **fuzzy** — `rapidfuzz token_sort_ratio ≥ FUZZY_THRESHOLD` (default 85), **gene-token
   guarded** to reject high-scoring wrong-gene siblings (e.g. "Signaling by BRCA1 mutants" vs
   "…AMER1 mutants");
4. otherwise → kept with a `NON-CANONICAL: ` prefix.

`scripts/build_synonyms.py` grows the synonym map offline (local-first; LLM only for genuinely
ambiguous names), optionally triggered after each run via `AUTO_UPDATE_SYNONYMS`.

### 3.6 Enrich (`run_enrich_stage`) — all non-LLM
- **STRING**: high-confidence partners (`combined_score ≥ STRING_MIN_SCORE`, default 700),
  attached as `string_interactors`.
- **GTEx safety** (two-tier): tier 1 sets `tier1_flag` when median TPM exceeds
  `GTEX_VITAL_TPM_THRESHOLD` (default 5) in **any** vital organ (brain, heart, liver, kidney,
  lung, adrenal); tier 2 sets `tier2_flag` when TPM exceeds `GTEX_TPM_THRESHOLD` (default 10)
  in ≥ `GTEX_TIER2_MIN_TISSUES` (default 2) secondary tissues. Records `safety_flag`,
  `tier1_flag`, `tier2_flag`, `tier1_high_tissues`/`tier2_high_tissues`, `max_vital_tpm`,
  `max_tpm`, the per-tier `safety_penalty` (0.60 / 0.80 / 1.0), and a review note.
- **CellxGene** (`ENABLE_CELLXGENE`): mean expression per cell type in `CENSUS_TISSUE` (cell
  types below `CENSUS_MIN_CELLS` dropped); attaches `cellxgene_expression` and unions the top 5
  cell types into `cellular_states` with a `CellxGene: ` prefix.

### 3.7 Network & scoring (`src/network.py`)
A `MultiDiGraph` (so typed edges coexist) with target nodes carrying the full annotation, plus
STRING satellite interactor nodes. Edge types: `pathway_comembership` (≥1 shared pathway,
weight = #shared), `direct_interaction` (gene → named LLM interactor that is also a node), and
`string_interaction` (weighted by combined score). Centrality is computed on an **undirected,
simple projection** `nx.Graph(G)` so direction and parallel edges don't distort it; only target
nodes are scored.

---

## 4. The composite score

```
composite = ( w_betweenness · betweenness
            + w_degree      · degree
            + w_disease     · min(disease_score, 1.0)
            + w_druggability· druggability_bonus
            + w_cellxgene   · cellxgene_score )
          × confidence × safety_penalty
```

| Term | How it's computed | Default weight |
|---|---|---|
| `betweenness` | normalized betweenness centrality on the undirected projection | `0.25` |
| `degree` | degree centrality (undirected projection) | `0.15` |
| `disease_score` | Σ over disease associations whose name matches a scoring term: strong 1.0 / moderate 0.5 / weak 0.2. **Capped at 1.0 in the composite**, reported uncapped in the TSV | `0.35` |
| `druggability_bonus` | `0.2` if `druggability_notes` non-empty, else `0.0` | `0.10` |
| `cellxgene_score` | `1.0` if ≥3 measured cell types, `0.5` if ≥1, else `0.0` | `0.15` |

The five weights live in `config` and are **validated to sum to 1.0 at startup**.
`confidence` is the merged record's confidence (0–1). `safety_penalty` is the two-tier GTEx
multiplier — `GTEX_TIER1_PENALTY` (0.60) for a tier-1 vital-organ flag, `GTEX_TIER2_PENALTY`
(0.80) for a tier-2-only flag, else 1.0 — a deprioritization, never elimination.

---

## 5. Two-layer cache

Two independent layers under `CACHE_DIR` (`outputs/cache/`, shared across runs — not inside
the per-run directory), so curating pathway canonicalization refreshes output **without
re-fetching/re-extracting**:

| Layer | Path | Key digests | Invalidated by |
|---|---|---|---|
| **Raw** | `raw/{gene}_{extract_key}.json` | gene, disease context, extraction model, PubMed depth (`PUBMED_MAX_RESULTS` + `PUBMED_EXTRACT_LIMIT`), `EXTRACTION_MAX_WORDS`, `EXTRACT_PROMPT_VERSION` | gene / source / disease / extraction-model / PubMed-depth / word-budget / prompt-version change |
| **Final** | `final/{gene}_{full_key}.json` | `extract_key` + synonym-file hash + Reactome-file hash + confidence threshold + merge model + enrich params | any of those, i.e. **also** synonym/reference/threshold/merge/enrich changes |

`make_cache_key` is a backward-compatible alias of `make_full_key`. File fingerprints are
memoized per `(path, mtime, size)` so the large Reactome file isn't re-hashed per gene.

**Execution order** (`run_gene`): final hit → raw hit (replay merge+enrich) → full chain.
A corrupt/truncated cache file is treated as a **miss** (logged WARNING, recomputed) — it can
never wedge a gene.

**Flags**: `ENABLE_CACHE` (off → no cache), `FORCE_RERUN` (bypass both, recompute all),
`FORCE_REMERGE` (bypass final only, replay merge+enrich from raw).

**Pruning** (`prune_stale_cache`, `PRUNE_CACHE` default on): after a successful run, for each
gene in the list it keeps only the current-key file in each layer and deletes the rest. Runs
**before** `AUTO_UPDATE_SYNONYMS` so a synonym rewrite can't retroactively mark the run's own
fresh files as stale.

---

## 6. Run accounting

Outcome counters are **disjoint**: `Total = Succeeded + Cached + Remerged + Failed`.
- `genes_succeeded` — ran the **full** fetch+extract+merge+enrich chain.
- `genes_cached` — final-cache hit (no work).
- `genes_remerged` — raw-cache hit (merge+enrich replayed, no fetch/extract).
- `genes_failed` — raised a retryable error past all retries (error-isolated, excluded).

The boxed run report also tracks LLM calls, token usage (incl. prompt-cache reads), estimated
cost, the `NON-CANONICAL` ratio, the `Cache pruned: X raw, Y final` line (only when > 0), and
runtime. Concurrency uses `asyncio.Semaphore`; counter increments are safe because asyncio is
single-threaded and no `await` splits an increment.

---

## 7. Configuration surface

All read once into `src/config.py` (`PipelineConfig`, re-reads env on each instantiation).
Grouped: models (`EXTRACTION_MODEL`, `MERGE_MODEL`, `MAX_TOKENS`, `EXTRACTION_MAX_WORDS`);
pipeline (`CONFIDENCE_THRESHOLD`, `PUBMED_MAX_RESULTS`, `PUBMED_EXTRACT_LIMIT`,
`SEMAPHORE_LIMIT`, `LOG_LEVEL`); disease (`DISEASE_CONTEXT`, `DISEASE_TERMS`); STRING
thresholds; two-tier GTEx (`GTEX_VITAL_TPM_THRESHOLD`, `GTEX_TPM_THRESHOLD`,
`GTEX_TIER2_MIN_TISSUES`, `GTEX_TIER1_PENALTY`, `GTEX_TIER2_PENALTY`); CellxGene
(`ENABLE_CELLXGENE`, `CENSUS_*`); scoring weights; fuzzy/synonym (`FUZZY_THRESHOLD`,
`SYNONYM_*`, `AUTO_UPDATE_SYNONYMS`); cache (`ENABLE_CACHE`, `CACHE_DIR`, `PRUNE_CACHE`,
`FORCE_RERUN`, `FORCE_REMERGE`); plot layout. See the README env table for defaults.

---

## 8. Caveats

These shape how outputs should be interpreted. Several were confirmed by live end-to-end runs.

1. **LLM non-determinism (reduced, not eliminated).** Extraction and merge are generative but
   both run at **`temperature=0`** (greedy), so most run-to-run drift is gone. temp-0 is
   *near*-deterministic, not bit-identical, so pathway sets, the `NON-CANONICAL` count, and
   rarely the exact ranking can still shift between reruns on the same genes. In normal
   operation the two-layer cache also freezes results — variance only appears when a stage
   actually re-runs (a `FORCE_RERUN`/`EXTRACT_PROMPT_VERSION` bump for extraction; any
   final-cache-invalidating change for merge). Treat any single run as a (now much more stable)
   sample, not ground truth.
2. **"Re-merge" is not zero-cost for multi-source genes.** Canonicalization is local, but it
   runs *after* the merge LLM call, and merge calls the model for any gene with ≥2 sources (the
   common case). So a synonym-edit rerun re-bills merge tokens; only single-source genes
   re-canonicalize for free. It is still far cheaper than a full rerun (no fetch/extract).
3. **`AUTO_UPDATE_SYNONYMS` makes the cache take two runs to settle.** A run can rewrite the
   synonym map after writing its finals, which shifts `full_key`; the *next* run then re-merges
   and only the run after that hits the final cache cleanly. Pruning is sequenced before the
   rewrite to avoid deleting the current run's own finals.
4. **Final cache is checked before raw.** Corrupting a raw file while the final entry is warm
   has no effect — the final hit short-circuits. Exercising raw-layer behavior requires the
   final entry to miss first.
5. **Pruning only covers genes in the current list.** Remove a gene from
   `inputs/target_genes.txt` and its cache files are never pruned (they're orphaned, not
   scanned). Stale files for *listed* genes are pruned every run.
6. **CellxGene census cache ignores `CENSUS_VERSION`.** The per-gene census file is keyed by
   `{gene}_{tissue}` only. Changing `CENSUS_VERSION` invalidates the *final* cache (it's in
   `full_key`) but the enrich step reuses the **stale** census data on disk. Delete
   `refs/census_cache/` to truly re-pull a new version.
7. **CellxGene first run is slow/expensive.** An uncached tissue query scans millions of cells
   (~10–20 min/gene). Results cache per `(gene, tissue)`; set `ENABLE_CELLXGENE=false` to skip.
8. **Network centrality needs scale (≥ ~20 genes).** A 5-gene graph is too sparse for
   betweenness/degree to carry signal; the network terms of the composite are noisy at small N.
9. **`EXTRACT_PROMPT_VERSION` is a manual lever.** It is *not* derived from the prompt text — if
   you change the extractor prompt (or decoding, e.g. temperature) without bumping the constant,
   stale raw extractions are reused. Currently `"2"` (v2 = extraction pinned to `temperature=0`;
   the bump invalidated v1 raw caches so extractions regenerate under greedy decoding).
10. **GTEx flags are a prompt for review, not a verdict.** High normal-tissue expression can be
    on-target biology (e.g. a tissue-specific TF like FOXF1 in lung). If the GTEx table can't
    be downloaded, the filter degrades to "no concern" (no flags) rather than failing.
11. **`disease_score` is uncapped in the TSV** but capped at 1.0 inside the composite — the two
    columns will disagree for highly-studied genes (e.g. TP53).
12. **`source_pmids` reflect inputs, not the merge model** (deliberate, to avoid fabricated
    citations) — they are the union of validated per-source PMIDs.
13. **Symbol resolution is HGNC-centric.** Non-canonical aliases may fail OpenTargets/STRING
    resolution and silently reduce a gene to fewer sources.
14. **Skipped genes fall outside the outcome buckets.** A gene with no source clearing the
    confidence gate returns `None` — it's neither succeeded, cached, remerged, nor failed, so
    `Total` can exceed the sum when skips occur.

---

## 9. Known gaps & differences

- **Batch mode parity.** `batch_pipeline.py` shares the cached extraction prompt but does
  **not** attach STRING/GTEx/CellxGene, so batch-built networks lack those edges, the safety
  penalty, and the `cellxgene_score` term.
- **No cross-config GC.** Pruning is per-listed-gene; there is no global sweep for files of
  genes no longer in the list.
- **Cytoscape export is a script**, not a `network.py` function — `scripts/export_cytoscape.py`
  emits Cytoscape.js JSON + CX2 (full + targets-only).

---

## 10. Where to look in the code

| Concern | Location |
|---|---|
| Orchestration, cache, counters, report | `pipeline.py` |
| Config (all env vars, weight validation) | `src/config.py` |
| LLM call (forced tool use + prompt caching) | `src/llm.py` |
| Pathway exact/synonym/fuzzy + gene-token guard | `src/pathways.py` |
| Per-source extraction + tool schema | `src/extractor.py` |
| Multi-source merge + canonicalization | `src/merger.py` |
| Graph build + composite scoring | `src/network.py` |
| Fetchers | `src/fetchers/{pubmed,uniprot,opentargets,string_db,cellxgene}.py` |
| GTEx safety | `src/filters/gtex_safety.py` |
| Synonym builder / Cytoscape export | `scripts/` |
| Tests (68) | `tests/` |
