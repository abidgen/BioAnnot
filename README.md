# Biological Annotation Agentic Pipeline

![Python](https://img.shields.io/badge/python-3.11-blue.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)

## Overview

This pipeline ingests a list of gene/protein targets and produces structured JSON
annotations — target functions, cellular states, pathway memberships, disease/biomarker
associations, interactors, and druggability notes — by querying PubMed, UniProt,
OpenTargets, and Reactome. It uses the Anthropic API with forced tool use to extract
structured annotations from each source and to merge them with conflict resolution, then
builds a NetworkX graph from the merged annotations to score and rank targets for
prioritization.

## Architecture

```
inputs/target_genes.txt
        │
        ▼
┌──────────────────────────────────────────────────┐
│  fetchers/  (async, per gene)                    │
│   ├─ pubmed.py       → abstracts + PMIDs         │
│   ├─ uniprot.py      → function, GO, locations   │
│   ├─ opentargets.py  → pathways, disease assoc.  │
│   ├─ reactome.py     → canonical pathway names   │
│   └─ string_db.py    → PPI partners (≥700)       │
└──────────────────────────────────────────────────┘
        │  raw source text (per gene, per source)
        │  [string_db output skips the LLM — see below]
        ▼
┌───────────────────────────────────────────────┐
│  extractor.py   (claude-opus-4-8, tool use)   │
│   one structured annotation per source        │
└───────────────────────────────────────────────┘
        │  list of per-source annotations
        ▼
┌───────────────────────────────────────────────┐
│  merger.py      (claude-sonnet-4-6, tool use) │
│   reconcile sources, resolve conflicts,       │
│   flag non-canonical pathways vs Reactome ref │
└───────────────────────────────────────────────┘
        │  one merged annotation per gene
        ▼
┌──────────────────────────────────────────────────────┐
│  filters/gtex_safety.py                              │
│   flag high normal-tissue expression (GTEx v8) →     │
│   attach safety_assessment to each merged record     │
└──────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────┐
│  network.py     (NetworkX)                           │
│   build graph (pathway_comembership, direct_inter-   │
│   action, string_interaction edges; + STRING         │
│   satellite nodes) + compute priority scores         │
│   (GTEx-flagged targets get a 0.75 composite penalty)│
└──────────────────────────────────────────────────────┘
        │
        ▼
   outputs/  (annotations.jsonl, final_annotations.json,
              target_network.gpickle, prioritized_targets.tsv)
```

The orchestrator [`pipeline.py`](pipeline.py) drives this flow with a concurrency limit of
3 genes at a time. All intermediate results are written to disk under `outputs/` — that
directory is the single source of truth between stages. STRING is the one fetcher whose
output is factual rather than free text, so its PPI partners **bypass the LLM extractor and
merger** and are attached to each merged record directly, feeding `network.py` as
`string_interaction` edges (and satellite interactor nodes). Likewise, the GTEx safety
filter ([`src/filters/gtex_safety.py`](src/filters/gtex_safety.py)) is a lookup, not an LLM
call: it attaches a `safety_assessment` to each merged record so `network.py` can
deprioritize targets that are highly expressed in sensitive normal tissues.

## Quick Start

From a clean checkout:

```bash
# 1. Create and activate the environment
mamba create -n bio_annot python=3.11 -y
mamba activate bio_annot

# 2. Install Python dependencies
pip install -r requirements.txt

# 3. Configure secrets (see below for what each key is and where to get it)
cp env.example .env
$EDITOR .env

# 4. Download the canonical Reactome pathway reference (human only)
mkdir -p refs
curl -s "https://reactome.org/download/current/ReactomePathways.txt" \
  | awk -F'\t' '$3=="Homo sapiens" {print $2}' \
  > refs/reactome_pathways.txt

# 5. Provide a gene list (one HGNC symbol per line)
mkdir -p inputs
printf "FOXF1\nTP53\nEGFR\nKRAS\nBRCA1\n" > inputs/target_genes.txt

# 6. Run the pipeline
python pipeline.py
```

### Required environment variables (`.env`)

Fill these into your `.env` file (never commit it — it is git-ignored):

| Key | Required | Purpose | Where to get it |
|---|---|---|---|
| `ANTHROPIC_API_KEY` | Yes | Authenticates the Anthropic API used by the extractor and merger | <https://console.anthropic.com/> → API Keys |
| `NCBI_EMAIL` | Yes | NCBI Entrez policy requires a contact email on every request | Your own email address |
| `NCBI_API_KEY` | Optional | Raises the NCBI rate limit from 3 → 10 requests/sec | NCBI account → Settings → API Key Management |
| `CONFIDENCE_THRESHOLD` | Optional | Drops extractions below this confidence (default `0.65`) | — |
| `DISEASE_CONTEXT` | Optional | Single disease label that focuses the run (default `cancer`; e.g. `fibrosis`, `neurodegeneration`) | — |
| `DISEASE_TERMS` | Optional | Comma-separated synonym list for the context (default `cancer,tumor,carcinoma,sarcoma,lymphoma,leukemia`) | — |
| `LOG_LEVEL` | Optional | Logging verbosity (default `INFO`) | — |

### Disease context

The run is no longer hardcoded to oncology. `DISEASE_CONTEXT` sets a single label and
`DISEASE_TERMS` a comma-separated synonym list; together they are resolved once by
`utils.load_disease_context()` and wired into three places: the PubMed search query (a short
OR-clause built from the context plus the first couple of terms and a generic `disease`
catch-all), the extractor's system prompt, and the prioritization scoring (a disease
association counts toward `disease_score` when any term matches its name). To retarget the
pipeline at, say, fibrosis, set `DISEASE_CONTEXT=fibrosis` and
`DISEASE_TERMS=fibrosis,fibrotic,scarring` — no code changes. This replaces the old
single-term `DISEASE_FILTER` variable.

## Repository Layout

```
bio-annotation-pipeline/
├── README.md                   ← this file
├── requirements.txt
├── env.example
│
├── inputs/
│   └── target_genes.txt        ← one gene symbol per line (e.g. FOXF1, TP53, EGFR)
│
├── refs/
│   ├── reactome_pathways.txt   ← canonical Reactome pathway names (one per line)
│   ├── uniprot_surface.txt     ← surface proteome gene list (optional filter)
│   └── gtex_median_tpm.gct.gz  ← GTEx v8 median-TPM table (auto-downloaded, cached)
│
├── src/
│   ├── fetchers/
│   │   ├── pubmed.py           ← PubMed/Entrez abstract fetcher
│   │   ├── uniprot.py          ← UniProt REST API fetcher
│   │   ├── opentargets.py      ← OpenTargets GraphQL fetcher
│   │   ├── reactome.py         ← Reactome pathway fetcher
│   │   └── string_db.py        ← STRING PPI interaction-partner fetcher
│   ├── filters/
│   │   └── gtex_safety.py      ← GTEx normal-tissue expression safety filter
│   ├── extractor.py            ← Anthropic API tool-use extraction core
│   ├── merger.py               ← LLM-assisted multi-source merge & conflict resolution
│   ├── network.py              ← NetworkX graph builder + target prioritization scorer
│   └── utils.py                ← logging, retry decorator, PMID validator
│
├── pipeline.py                 ← main orchestrator (run this)
├── batch_pipeline.py           ← Anthropic Batch API variant for 50+ genes
├── visualize_network.py        ← plots from existing outputs (no rerun)
│
└── outputs/                    ← auto-created at runtime
    ├── raw/                    ← per-gene per-source raw extraction JSONs
    ├── annotations.jsonl       ← merged annotation per gene (newline-delimited JSON)
    ├── final_annotations.json  ← full merged dict keyed by gene symbol
    ├── target_network.gpickle  ← NetworkX graph
    ├── prioritized_targets.tsv ← ranked target table
    └── plots/                  ← visualization PNGs (after visualize_network.py)
```

## Output Files

All outputs land under `outputs/` and are regenerated on each run.

- **`outputs/raw/{gene}_raw.json`** — the per-source annotation records (PubMed, UniProt,
  OpenTargets) for one gene before merging. Useful for debugging where an annotation came
  from or why a source was dropped below the confidence threshold.

- **`outputs/annotations.jsonl`** — one merged annotation per line (newline-delimited
  JSON), appended as each gene completes. Convenient for streaming/`jq` processing and as
  an incremental record even if a later gene fails.

- **`outputs/final_annotations.json`** — the full merged dictionary keyed by gene symbol,
  written once at the end. This is the canonical structured result. Each value contains
  `functions`, `cellular_states`, `pathways` (non-canonical names prefixed
  `NON-CANONICAL: `), `disease_associations` (each with `role` and `evidence_strength`),
  `interactors`, `druggability_notes`, `confidence`, `source_count`, `source_pmids`, and
  `merged_at`.

- **`outputs/prioritized_targets.tsv`** — the ranked target table, one row per gene sorted
  by `composite` descending. Columns: `gene`, `composite`, `betweenness`, `degree`,
  `disease_score`, `druggability_bonus`, `confidence`, `safety_flag`,
  `safety_penalty_applied`, `high_expression_tissues`, `max_tpm`, `pathways`,
  `disease_associations` (list/dict fields flattened to pipe-separated strings). The
  `composite` score combines network centrality, disease relevance, and druggability, scaled
  by extraction confidence; targets flagged by the GTEx safety filter are then scaled by an
  additional `0.75` penalty (`safety_penalty_applied=True`), deprioritizing rather than
  eliminating them.

- **`outputs/target_network.gpickle`** — the NetworkX graph (a `MultiDiGraph`) pickled to
  disk. It holds the 5 (or however many) target nodes (`node_type="target"`, carrying the
  full annotation as attributes) plus satellite interactor nodes (`node_type="interactor"`)
  — the STRING partners of the targets, added so otherwise-isolated genes gain connectivity.
  Edge types are `pathway_comembership` (genes sharing ≥1 pathway), `direct_interaction`
  (a gene to a named LLM-extracted interactor that is also a node), and `string_interaction`
  (STRING PPI partners, weighted by `combined_score` on the 0–1000 scale — between two
  targets, or from a target to a satellite). Centrality is computed over the whole graph but
  only target nodes are scored and written to the TSV. Pass
  `build_target_network(..., include_interactor_nodes=False)` for a target-only graph with
  no satellites (cleaner for larger gene sets). Load with `pickle.load(open(path, "rb"))`.

## Quality Gates

These are enforced automatically by the pipeline:

| Gate | Rule | Action |
|---|---|---|
| Confidence filter | Drop extractions < 0.65 (`CONFIDENCE_THRESHOLD`) | Log warning, skip source |
| PMID validation | Only digits, 7–8 chars | Drop invalid, log |
| Pathway canonicity | Check against Reactome reference set | Prefix with `NON-CANONICAL: ` |
| Source agreement | Pathway needs ≥2 sources unless confidence ≥ 0.85 | Merger rule |
| Normal-tissue safety | >10 TPM in ≥3 sensitive GTEx tissues | Flag and apply 0.75 composite penalty (deprioritize) |
| Rate limiting | Max 3 concurrent gene fetches | `asyncio.Semaphore(3)` |

## Batch Mode

Use [`batch_pipeline.py`](batch_pipeline.py) instead of `pipeline.py` when processing
**≥ 50 genes**. It submits all extraction requests through the Anthropic Batch API, which
delivers roughly **50% cost reduction** in exchange for asynchronous (polled) completion
rather than real-time results. For small gene sets the standard `pipeline.py` is simpler
and returns faster.

```bash
python batch_pipeline.py
```

The batch job ID is written to `outputs/batch_id.txt`; the script polls until the batch
ends, then collects results and runs the same merge → network → output stages as the
standard pipeline.

## Visualization

After a pipeline run has produced `outputs/target_network.gpickle` and
`outputs/prioritized_targets.tsv`, generate plots with:

```bash
python visualize_network.py
```

This reads existing outputs only (no pipeline rerun) and writes three PNGs to
`outputs/plots/`:

- **`target_network.png`** — the target graph with nodes colored and sized by composite
  score (red = high, blue = low), edges colored by type (green = pathway co-membership,
  orange = direct interaction), edge type labels, and a composite-score colorbar.
- **`score_breakdown.png`** — per-gene horizontal bars comparing confidence,
  disease_score (capped at 1.0, its effective contribution to the composite),
  network_score (betweenness + degree), and composite_score on a shared 0–1 scale.
- **`pathway_heatmap.png`** — a genes × canonical-pathways presence matrix
  (1 = gene has pathway). Non-canonical pathway names are excluded.

## Extension Points

Two enrichment layers are now **implemented**:

- **STRING PPI enrichment** (`src/fetchers/string_db.py`) — see the Architecture section and
  the `string_interaction` edges in `outputs/target_network.gpickle`.
- **GTEx safety filter** (`src/filters/gtex_safety.py`) — flags genes highly expressed in
  sensitive normal tissues (>10 TPM in ≥3 tissues) and applies a 0.75 composite penalty in
  `network.py`. The GTEx v8 median-TPM table is auto-downloaded and cached at
  `refs/gtex_median_tpm.gct.gz` on first use.

The following extensions are still planned:

- **CellxGene Census** — `src/fetchers/cellxgene.py` fetching mean expression per cell type
  for each gene (e.g. in human lung tissue) to ground the `cellular_states` field.
- **Cytoscape export** — a `network.py` `export_cytoscape_json(G, path)` helper using
  `nx.cytoscape_data(G)` for interactive exploration.
- **Batch-mode parity** — `batch_pipeline.py` does not yet attach STRING partners or GTEx
  safety assessments to records, so batch-built networks currently lack those edges and
  penalties.

## Known Limitations

- **LLM non-determinism causes run-to-run variance.** The extractor and merger are
  generative, so repeated runs on the same genes can return different pathway sets and
  slightly different phrasing. Counts (e.g. number of `NON-CANONICAL` flags) and even the
  exact ranking can shift between runs; treat individual runs as samples, not fixed truth.

- **Remaining `NON-CANONICAL` flags are mostly informal signaling names.** After
  normalization (case-insensitive matching plus stripping Reactome `(R-HSA-…)` stable-ID
  suffixes), the pathways still flagged are typically informal shorthand the model emits —
  e.g. `PI3K/AKT/mTOR signaling`, `RAS-RAF-MEK-ERK (MAPK) cascade`, `JAK/STAT signaling` —
  rather than exact Reactome names. These are genuine non-canonical names, not a matching
  bug; the gate is working as intended.

- **Network centrality is only meaningful on larger gene sets (≥ ~20 genes).** With a
  handful of genes the graph is too sparse for betweenness/degree centrality to carry
  signal (a 5-gene network often has just a few edges), so the network-derived component of
  the composite score is noisy at small scale. Run with a substantial target list before
  relying on the centrality terms.

- **`source_pmids` are unioned from inputs, not taken from the merger model.** The
  `annotate_target` tool schema has no `source_pmids` field, so the merged record's PMIDs
  are computed by unioning the validated PMIDs from the per-source inputs rather than read
  back from the merge model's output. This is deliberate (it prevents the model from
  inventing citations) but means the merged `source_pmids` reflect the input sources, not
  any model-level attribution.

## License

Released under the [MIT License](LICENSE).

