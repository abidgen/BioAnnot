# Information Flow — From Gene List to Ranked Targets

> A top-to-bottom trace of how data moves through the pipeline, using the real
> 5-gene validation run in `outputs/` (`FOXF1, TP53, EGFR, KRAS, BRCA1`).
> Every number below is taken from the actual run logs and output files, not
> illustrative placeholders.

---

## 0. The input

`inputs/target_genes.txt` — one HGNC symbol per line:

```
FOXF1
TP53
EGFR
KRAS
BRCA1
```

`load_gene_list()` strips blanks/comments and uppercases. These 5 symbols are
processed concurrently (max 3 at a time, `asyncio.Semaphore(3)`).

The rest of this document follows **FOXF1** end-to-end, then shows how the
**ranking** across all five genes is computed.

---

## 1. Fetch — four public sources queried per gene

For each gene the pipeline fires off independent API calls. Here is exactly what
happened for FOXF1 (from `outputs/pipeline.log`):

| Source | Query sent | What came back |
|---|---|---|
| **PubMed** (Entrez) | `esearch` `term=FOXF1[gene] AND (cancer OR tumor OR disease)`, `retmax=20`, then `efetch` the abstracts | **20 PMIDs** → 20 abstracts (title + abstract + year + journal) |
| **UniProt** | `gene_exact:FOXF1 AND organism_id:9606 AND reviewed:true` | accession **Q12946**, 32 GO terms, 1 subcellular location, function text |
| **OpenTargets** | GraphQL: resolve symbol → `ENSG00000103241`, then fetch target associations | **2 pathways, 10 disease associations**, function descriptions |
| **STRING** | `interaction_partners?identifiers=FOXF1&species=9606&required_score=700` | **10 partners** (SHH 915, BMP4 886, GATA4 871, …) |

Two more data sources attach later (they don't go through the LLM):

- **GTEx** (`refs/gtex_median_tpm.gct.gz`) — normal-tissue expression for the safety filter.
- **CellxGene Census** — measured per-cell-type expression in lung tissue (170 cell types for FOXF1).

### PMID safety
PMIDs are passed through `validate_pmids()` (regex `^\d{7,8}$`). The log line
`search_pmids(FOXF1): 20 PMIDs after validation` confirms all 20 survived.
**The LLM never emits PMIDs** — they are attached programmatically downstream so
they cannot be hallucinated.

---

## 2. Extract — LLM turns free text into structured records

This is the core LLM step (`src/extractor.py` + `src/llm.py`). Each source's
text is sent **separately** to the model, producing one structured record per
source.

**What the LLM does:** it is given a single tool, `annotate_target`, and called
with `tool_choice={"type":"tool","name":"annotate_target"}` — this **forces** the
model to return a JSON object matching the schema (no free-text parsing). The
model fills in `functions`, `cellular_states`, `pathways`,
`disease_associations` (each with `role` + `evidence_strength`), `interactors`,
`druggability_notes`, and a `confidence` number.

- **Model:** `claude-opus-4-8` for extraction.
- **System prompt** (abridged): *"senior biomedical annotation scientist…
  extract only what the text states, do not import outside knowledge…
  use exact Reactome pathway names… never invent interactors or PMIDs."*
  The active disease context (`cancer`) is injected so extraction stays on-topic.
- **Truncation:** input text capped at ~3000 words to avoid context overflow.
- **Caching:** the static system block + tool schema are marked
  `cache_control: ephemeral`, so they're billed once and read from cache after.
  FOXF1 PubMed call logged: `input=6063 output=1104 cache_read=0 cache_created=1339`,
  then later calls show `cache_read=1339`.

### How confidence is measured (the key question)

`confidence` is **produced by the LLM itself**, not computed by code. The system
prompt gives the model an explicit rubric:

```
0.9–1.0  Multiple strong sources, consistent findings
0.7–0.9  Two or more sources, minor inconsistencies
0.5–0.7  Single source, or conflicting evidence
0.0–0.5  Weak/indirect evidence, or text only tangentially about the gene
```

The model judges *how well the supplied text supports a high-quality annotation
for this specific gene* and emits a single 0.0–1.0 number. It is told to prefer
**precision over recall** — return fewer items and a lower confidence rather than
padding a sparse record.

**FOXF1's three extracted records:**

| Source | confidence | functions | pathways | diseases |
|---|---|---|---|---|
| PubMed (20 rich abstracts) | **0.88** | 8 | 8 | 9 |
| UniProt (terse curated entry) | **0.55** | 7 | 2 | 1 |
| OpenTargets (associations only) | **0.55** | 2 | 2 | 10 |

The model rated the dense PubMed text high (0.88) and the thin UniProt/OpenTargets
text low (0.55) — exactly as the rubric intends.

Raw per-source records are saved to `outputs/raw/FOXF1_raw.json`.

---

## 3. Filter — the confidence gate

`pipeline.py` keeps only sources at or above the threshold
(`CONFIDENCE_THRESHOLD = 0.65`):

```python
sources = [a for a in [pubmed_ann, uniprot_ann, ot_ann]
           if a and a.get("confidence", 0) >= config.confidence_threshold]
```

For **FOXF1**: PubMed (0.88) passes; UniProt (0.55) and OpenTargets (0.55) are
**dropped**. → **only 1 source survives.** This is why FOXF1's final record shows
`"source_count": 1`. If *zero* sources passed, the gene would be logged
(`No high-confidence sources for …`) and skipped entirely.

---

## 4. Merge — reconcile surviving sources into one record

`src/merger.py`, model **`claude-sonnet-4-6`** (cheaper; merge is easier than
extraction), again with the forced `annotate_target` tool.

- **Single-source genes (like FOXF1):** no reconciliation needed — the merger
  short-circuits, returns the record with `source_count=1`, and just validates
  pathway names.
- **Multi-source genes:** the LLM applies the merge rules:
  - Include a pathway only if **≥2 sources agree** OR a single source rates
    confidence **≥ 0.85**.
  - Conflicting disease roles (oncogene vs tumor_suppressor) → note
    context-dependence, don't pick one.
  - Union interactors across sources.
  - **Merged confidence** = reliability-weighted average of source confidences
    (curated UniProt/OpenTargets weigh more than free-text PubMed), minus 0.1 per
    material conflict resolved, clamped to [0, 1].

### Pathway canonicity
After the merge, every pathway is checked against `refs/reactome_pathways.txt`.
Non-matching names are **prefixed** (not dropped) with `NON-CANONICAL: `. In
FOXF1's final record, `JAK-STAT signaling` became
`NON-CANONICAL: JAK-STAT signaling`, while `Signaling by WNT`, `Cellular
Senescence`, etc. matched and were kept clean.

### Non-LLM enrichment attached here
- **STRING** 10 partners → `string_interactors`.
- **GTEx safety** → `safety_assessment`. FOXF1 flagged: high expression in
  Colon-Sigmoid, Lung, Small Intestine (max **122.1 TPM**). Logged as a WARNING.
- **CellxGene** → `cellxgene_expression` (top cell types: pericyte, endothelial,
  …; 170 cell types total).

Result is written to `outputs/final_annotations.json` (and appended to
`annotations.jsonl`). FOXF1's merged record: 8 functions, 8 pathways
(1 non-canonical), 9 disease associations, 9 interactors unioned with 10 STRING
partners, `confidence: 0.88`.

---

## 5. Network — build the graph

`src/network.py` builds a `MultiDiGraph` over the 5 genes (+ STRING satellite
nodes):

- **Nodes:** one per gene, carrying all merged attributes.
- **Edges:**
  - `pathway_comembership` between genes sharing ≥1 pathway (weight = #shared).
  - `direct_interaction` for named interactors that are also nodes.
  - `string_interaction` to STRING partners (weighted by combined score/1000).

Centrality (betweenness, degree) is computed on the **undirected projection**
`nx.Graph(G)` so edge direction/multiplicity don't distort it. Saved to
`outputs/target_network.gpickle`.

---

## 6. Score & rank — the composite

For every gene (`compute_priority_scores`):

```
disease_score = Σ  { strong:1.0, moderate:0.5, weak:0.2 }
                   over disease_associations whose name matches a disease term

composite = ( 0.25·betweenness
            + 0.15·degree
            + 0.35·min(disease_score, 1.0)
            + 0.10·druggability_bonus      (0.2 if any druggability_notes)
            + 0.15·cellxgene_score )       (1.0 ≥3 cell types, 0.5 ≥1, else 0)
            × confidence
            × safety_penalty               (0.75 if GTEx-flagged, else 1.0)
```

### Worked example — why BRCA1 beats TP53

| | betweenness | degree | disease (capped) | drug | cellx | conf | safety | **composite** |
|---|---|---|---|---|---|---|---|---|
| **BRCA1** | 0.397 | 0.268 | 1.0 | 0.2 | 1.0 | 0.88 | ×1.00 | **0.580** |
| **TP53** | 0.546 | 0.274 | 1.0 | 0.2 | 1.0 | 0.89 | ×0.75 | **0.466** |

TP53 has **higher** centrality and confidence, yet ranks **below** BRCA1. The
reason is the **GTEx safety penalty**: TP53 is broadly expressed in normal
tissues (flagged), so its composite is multiplied by 0.75; BRCA1 is not flagged
(×1.0). Check: BRCA1 = (0.25·0.397 + 0.15·0.268 + 0.35 + 0.02 + 0.15) × 0.88 =
**0.5806** ✓. The disease_score raw value (BRCA1 7.2, TP53 14.4) is **capped at
1.0** in the composite, so being "even more cancer-associated" gives TP53 no
extra credit.

### Final ranking (`outputs/prioritized_targets.tsv`)

```
1. BRCA1  0.580   (not safety-flagged)
2. TP53   0.466   (safety-flagged ×0.75)
3. EGFR   0.459   (safety-flagged ×0.75)
4. KRAS   0.428   (safety-flagged ×0.75)
5. FOXF1  0.380   (lowest centrality: degree 0.068, only 1 source)
```

FOXF1 lands last: as a tissue-specific transcription factor it has the fewest
graph connections (degree 0.068) and survived on a single source, so its
centrality contribution is small — even though its individual annotation is rich.

---

## Flow-of-information diagram

```
 inputs/target_genes.txt   [FOXF1 TP53 EGFR KRAS BRCA1]
            │
            ▼   load_gene_list()  → uppercase, ≤3 concurrent
 ┌──────────────────────── per gene (e.g. FOXF1) ──────────────────────────┐
 │                                                                         │
 │   1. FETCH (parallel API calls)                                         │
 │   ┌─────────┬──────────┬────────────┬─────────┐                         │
 │   │ PubMed  │ UniProt  │ OpenTargets│ STRING  │   + GTEx  + CellxGene   │
 │   │ 20 PMIDs│ Q12946   │ ENSG…03241 │ 10 ptnrs│   (non-LLM enrichment)  │
 │   │ abstrac.│ 32 GO    │ 2 pw/10 dis│         │                         │
 │   └────┬────┴────┬─────┴─────┬──────┴────┬────┘                         │
 │        │         │           │           │                              │
 │        ▼ validate_pmids()    │           │                              │
 │   2. EXTRACT  (claude-opus-4-8, forced annotate_target tool)            │
 │        │         │           │           │                              │
 │     conf 0.88  conf 0.55   conf 0.55     │ ← LLM self-rates per rubric  │
 │        │         │           │           │                              │
 │        ▼─────────▼───────────▼           │   outputs/raw/FOXF1_raw.json │
 │   3. FILTER  (keep conf ≥ 0.65)          │                              │
 │     KEEP ✓    DROP ✗      DROP ✗         │                              │
 │        │                                 │                              │
 │        ▼                                 │                              │
 │   4. MERGE  (claude-sonnet-4-6)  ◄───── STRING / GTEx / CellxGene join  │
 │     source_count=1 · pathway canonicity check (NON-CANONICAL: …)        │
 │        │                                                                │
 │        ▼   outputs/final_annotations.json                               │
 └────────┼────────────────────────────────────────────────────────────────┘
          │  (all 5 genes merged)
                ▼
   5. NETWORK   build MultiDiGraph: pathway_comembership / direct /
                string_interaction edges → target_network.gpickle
          │
                ▼
   6. SCORE     composite = (0.25·btwn + 0.15·deg + 0.35·disease
                + 0.10·drug + 0.15·cellx) × confidence × safety_penalty
          │
                ▼
   outputs/prioritized_targets.tsv
   ┌─────────────────────────────────────┐
   │ 1 BRCA1 0.580   (no penalty)        │
   │ 2 TP53  0.466   (×0.75 safety)      │
   │ 3 EGFR  0.459   (×0.75 safety)      │
   │ 4 KRAS  0.428   (×0.75 safety)      │
   │ 5 FOXF1 0.380   (low degree, 1 src) │
   └─────────────────────────────────────┘
```

---

## One-line summary per stage

| Stage | Input | Output | Key filter / rule |
|---|---|---|---|
| Fetch | gene symbol | PMIDs, UniProt JSON, OT associations, STRING partners | PMID regex `^\d{7,8}$` |
| Extract | source text | structured record + **LLM-rated confidence** | force tool, ~3000-word cap |
| Filter | 3 records | high-conf subset | `confidence ≥ 0.65` |
| Merge | surviving records | one merged record | ≥2 sources or conf ≥0.85 per pathway; NON-CANONICAL prefix |
| Network | all merged records | graph | centrality on undirected projection |
| Score | graph | ranked TSV | disease_score capped at 1.0; safety ×0.75 |
```
