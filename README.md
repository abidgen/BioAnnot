# Biological Annotation Agentic Pipeline

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
│   └─ reactome.py     → canonical pathway names   │
└──────────────────────────────────────────────────┘
        │  raw source text (per gene, per source)
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
┌─────────────────────────────────────────────────┐
│  network.py     (NetworkX)                      │
│   build graph + compute priority scores         │
└─────────────────────────────────────────────────┘
        │
        ▼
   outputs/  (annotations.jsonl, final_annotations.json,
              target_network.gpickle, prioritized_targets.tsv)
```

The orchestrator [`pipeline.py`](pipeline.py) drives this flow with a concurrency limit of
3 genes at a time. All intermediate results are written to disk under `outputs/` — that
directory is the single source of truth between stages.

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
| `DISEASE_FILTER` | Optional | Disease term used for prioritization scoring (default `cancer`) | — |
| `LOG_LEVEL` | Optional | Logging verbosity (default `INFO`) | — |

## Repository Layout

```
bio-annotation-pipeline/
├── README.md                   ← this file
├── requirements.txt
├── .env.example
│
├── inputs/
│   └── target_genes.txt        ← one gene symbol per line (e.g. FOXF1, TP53, EGFR)
│
├── refs/
│   ├── reactome_pathways.txt   ← canonical Reactome pathway names (one per line)
│   └── uniprot_surface.txt     ← surface proteome gene list (optional filter)
│
├── src/
│   ├── fetchers/
│   │   ├── pubmed.py           ← PubMed/Entrez abstract fetcher
│   │   ├── uniprot.py          ← UniProt REST API fetcher
│   │   ├── opentargets.py      ← OpenTargets GraphQL fetcher
│   │   └── reactome.py         ← Reactome pathway fetcher
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
  `disease_score`, `druggability_bonus`, `confidence`, `pathways`, `disease_associations`
  (list/dict fields flattened to pipe-separated strings). The `composite` score combines
  network centrality, disease relevance, and druggability, scaled by extraction
  confidence.

- **`outputs/target_network.gpickle`** — the NetworkX graph (a `MultiDiGraph`) pickled to
  disk. Nodes are genes carrying the full annotation as attributes; edges are
  `pathway_comembership` (genes sharing ≥1 pathway) and `direct_interaction` (a gene to a
  named interactor that is also a node). Load with
  `pickle.load(open(path, "rb"))`.

## Quality Gates

These are enforced automatically by the pipeline:

| Gate | Rule | Action |
|---|---|---|
| Confidence filter | Drop extractions < 0.65 (`CONFIDENCE_THRESHOLD`) | Log warning, skip source |
| PMID validation | Only digits, 7–8 chars | Drop invalid, log |
| Pathway canonicity | Check against Reactome reference set | Prefix with `NON-CANONICAL: ` |
| Source agreement | Pathway needs ≥2 sources unless confidence ≥ 0.85 | Merger rule |
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

The following extensions are planned:

- **STRING PPI** — `src/fetchers/string_db.py` querying STRING for interactor gene symbols
  with combined scores ≥ 700, to enrich the direct-interaction edges.
- **CellxGene Census** — `src/fetchers/cellxgene.py` fetching mean expression per cell type
  for each gene (e.g. in human lung tissue) to ground the `cellular_states` field.
- **GTEx safety filter** — `src/filters/gtex_safety.py` flagging genes with high normal
  tissue expression (>10 TPM in ≥3 sensitive tissues) as potential safety concerns.
- **Cytoscape export** — a `network.py` `export_cytoscape_json(G, path)` helper using
  `nx.cytoscape_data(G)` for interactive exploration.

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
```

