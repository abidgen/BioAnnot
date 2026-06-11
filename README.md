# Biological Annotation Agentic Pipeline

LLM-assisted structured extraction of biological annotations (target functions, cellular
states, pathway memberships, biomarker associations) from literature and public knowledge
bases, supporting target prioritization and network analysis.

See [CLAUDE.md](CLAUDE.md) for the full specification, module-by-module build steps, and
run instructions.

## Quick start

```bash
mamba create -n bio_annot python=3.11 -y && mamba activate bio_annot
pip install -r requirements.txt
cp .env.example .env   # fill in ANTHROPIC_API_KEY and NCBI_EMAIL
python pipeline.py
```

## Layout

```
src/fetchers/   PubMed, UniProt, OpenTargets, Reactome fetchers
src/extractor.py  Anthropic tool-use extraction core
src/merger.py     multi-source merge & conflict resolution
src/network.py    NetworkX graph builder + target prioritizer
pipeline.py       main orchestrator (< 50 genes)
batch_pipeline.py Batch API variant (>= 50 genes)
outputs/          auto-created results (raw JSONs, merged annotations, network, TSV)
```
