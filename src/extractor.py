"""Core LLM extraction module — Anthropic tool-use forced structured output.

Rules (per CLAUDE.md Step 5):
- Always use the extraction model and force the annotate_target tool.
- Never hallucinate PMIDs — only pass through IDs the fetcher returned.
- Truncate input text to ~3000 words before sending.
- Log token usage for every API call.
"""

from __future__ import annotations

import json
import logging
from typing import Any

from src.config import config
from src.llm import call_tool
from src.utils import (
    CACHE_MIN_TOKENS,
    estimate_tokens,
    load_disease_context,
    validate_pmids,
)

log = logging.getLogger("bio_annot.extractor")

# Model and token budget are centralized in src.config (env-configurable).
EXTRACTION_MODEL = config.extraction_model

MAX_INPUT_WORDS = 3000
MAX_TOKENS = config.max_tokens

# Detailed, static base prompt. Kept long and stable so that, together with the
# tool schema, it forms a cacheable prefix above the prompt-cache minimum. Only
# the small disease-context suffix (build_system_prompt) varies per run.
SYSTEM_PROMPT = (
    "You are a senior biomedical annotation scientist at the National Cancer "
    "Institute. Extract precise, evidence-based annotations ONLY from the "
    "provided text. Use canonical gene symbols (HGNC). Never assert anything the "
    "text does not support, and never speculate beyond it.\n\n"
    "You will receive text about a single gene/protein and must populate the "
    "annotate_target tool. Follow these field-by-field instructions exactly.\n\n"
    "FIELDS:\n"
    "- gene_symbol: The official HGNC symbol for the target, uppercase. Do not "
    "substitute an alias when the canonical symbol is known.\n"
    "- functions: Molecular and cellular functions (max 8 items). Be specific — "
    "prefer 'sequence-specific DNA-binding transcription factor' over 'regulates "
    "transcription'. Each item is a short phrase, not a sentence.\n"
    "- cellular_states: Cell types, tissues, or physiological/disease states in "
    "which the target is expressed or active. Use standard cell-type nomenclature "
    "where possible.\n"
    "- pathways: Signaling or metabolic pathways the target participates in. "
    "Always use exact Reactome pathway names. Valid examples: Signaling by WNT, "
    "RAF/MAP kinase cascade, Transcriptional Regulation by TP53. Do not invent "
    "pathway names; if the text mentions only an informal pathway, map it to the "
    "closest canonical Reactome name.\n"
    "- disease_associations: Each entry has a disease name, a role (oncogene, "
    "tumor_suppressor, biomarker, therapeutic_target, or unknown), and an "
    "evidence_strength (strong, moderate, weak). Assign a role only when the text "
    "supports it; otherwise use unknown.\n"
    "- interactors: Direct protein interactors EXPLICITLY named in the text, gene "
    "symbols only. Never invent interactors that are not named in the text.\n"
    "- druggability_notes: Note specific binding pockets, approved drugs, clinical "
    "trial agents, and resistance mechanisms mentioned in the text. Name the drug "
    "or drug class where given. Leave empty if the text says nothing about "
    "druggability.\n"
    "- confidence: A single number from 0.0 to 1.0 reflecting how well the text "
    "supports a high-quality annotation for THIS gene.\n\n"
    "PMID RULES:\n"
    "Never invent PMIDs. Only use IDs provided in the input. Source PMIDs are "
    "attached programmatically downstream — do not emit them yourself.\n\n"
    "CONFIDENCE RUBRIC:\n"
    "- 0.9-1.0: Multiple strong sources, consistent findings.\n"
    "- 0.7-0.9: Two or more sources, minor inconsistencies.\n"
    "- 0.5-0.7: Single source, or conflicting evidence.\n"
    "- 0.0-0.5: Weak or indirect evidence only, or text only tangentially related "
    "to the gene.\n\n"
    "GENERAL GUIDANCE:\n"
    "Extract only what the text states; do not import outside knowledge about the "
    "gene, however well known. Keep list entries concise and non-redundant, and "
    "deduplicate near-identical items. When the text is sparse or only tangentially "
    "about the gene, return fewer items and a lower confidence rather than padding "
    "the record. Prefer precision over recall: a short, accurate annotation is more "
    "useful downstream than a long, speculative one.\n\n"
    "Use exact Reactome pathway names wherever possible. Never invent interactors "
    "or pathways not mentioned in the text."
)


def build_system_prompt() -> str:
    """System prompt with the active disease context injected.

    Reads DISEASE_CONTEXT so extraction focuses on the configured disease area
    (cancer, fibrosis, neurodegeneration, …) rather than being hardcoded.
    """
    context = load_disease_context()["context"]
    return (
        f"{SYSTEM_PROMPT} "
        f"Focus annotations on relevance to {context}. "
        f"Prioritize disease associations, pathways, and cellular states "
        f"relevant to {context}."
    )

ANNOTATION_TOOL = [
    {
        "name": "annotate_target",
        "description": (
            "Extract structured biological annotations for a gene/protein from text"
        ),
        # Cache the (static) tool definition so it is not re-billed on every call.
        # This is the last/only tool, so the breakpoint covers the whole tools block.
        "cache_control": {"type": "ephemeral"},
        "input_schema": {
            "type": "object",
            "properties": {
                "gene_symbol": {"type": "string"},
                "functions": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Molecular and cellular functions, max 8 items",
                },
                "cellular_states": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Cell types or states where target is active/expressed",
                },
                "pathways": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Signaling or metabolic pathways — use Reactome/KEGG names",
                },
                "disease_associations": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "disease": {"type": "string"},
                            "role": {
                                "type": "string",
                                "enum": [
                                    "oncogene",
                                    "tumor_suppressor",
                                    "biomarker",
                                    "therapeutic_target",
                                    "unknown",
                                ],
                            },
                            "evidence_strength": {
                                "type": "string",
                                "enum": ["strong", "moderate", "weak"],
                            },
                        },
                        "required": ["disease", "role", "evidence_strength"],
                    },
                },
                "interactors": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Direct protein interactors mentioned, gene symbols only",
                },
                "druggability_notes": {
                    "type": "string",
                    "description": "Any mentions of druggability, binding pockets, drug classes",
                },
                "confidence": {
                    "type": "number",
                    "description": "0.0–1.0 extraction confidence given text quality and relevance",
                },
            },
            "required": ["gene_symbol", "functions", "pathways", "confidence"],
        },
    }
]

def _log_cache_prefix_size() -> None:
    """Confirm the static system+tools prefix is large enough to cache.

    The prompt cache only engages when the cached prefix clears CACHE_MIN_TOKENS,
    so we estimate it once at import. A short prefix is logged as a WARNING (which
    surfaces even before setup_logging configures handlers) so an undersized
    prompt is never silently un-cached.
    """
    system = build_system_prompt()
    prefix = system + json.dumps(ANNOTATION_TOOL)
    sys_tokens = estimate_tokens(system)
    prefix_tokens = estimate_tokens(prefix)
    level = logging.INFO if prefix_tokens >= CACHE_MIN_TOKENS else logging.WARNING
    log.log(
        level,
        "extractor system prompt ~%d tokens; system+tools cache prefix ~%d tokens "
        "(cache min %d)",
        sys_tokens,
        prefix_tokens,
        CACHE_MIN_TOKENS,
    )


_log_cache_prefix_size()


def _truncate_words(text: str, max_words: int = MAX_INPUT_WORDS) -> str:
    """Truncate text to at most max_words whitespace-delimited words."""
    words = text.split()
    if len(words) <= max_words:
        return text
    log.info("Truncating input from %d to %d words", len(words), max_words)
    return " ".join(words[:max_words])


async def extract_from_text(
    gene: str, text: str, source_pmids: list[str]
) -> dict[str, Any]:
    """Extract structured annotations for a gene from free text.

    Calls Claude via :func:`src.llm.call_tool` with the annotate_target tool
    forced (prompt caching + token logging handled there), then attaches
    validated source PMIDs to the returned annotation.
    """
    truncated = _truncate_words(text)

    user_content = f"Gene: {gene}\n\n{truncated}"

    annotation_input, _usage = await call_tool(
        model=EXTRACTION_MODEL,
        system_prompt=build_system_prompt(),
        messages=[{"role": "user", "content": user_content}],
        tools=ANNOTATION_TOOL,
        tool_name="annotate_target",
        max_tokens=MAX_TOKENS,
        label=f"extract_from_text({gene})",
    )

    annotation = dict(annotation_input)
    # Only pass through validated PMIDs — never invent or forward unvalidated IDs.
    annotation["source_pmids"] = validate_pmids(source_pmids)
    return annotation


def _format_uniprot_text(uniprot_data: dict[str, Any]) -> str:
    """Render a UniProt record dict as a readable text block."""
    lines = [
        f"Accession: {uniprot_data.get('accession', '')}",
        f"Gene name: {uniprot_data.get('gene_name', '')}",
        f"Protein name: {uniprot_data.get('protein_name', '')}",
        f"Function: {uniprot_data.get('function_text', '')}",
    ]
    locations = uniprot_data.get("subcellular_locations", [])
    if locations:
        lines.append("Subcellular locations: " + ", ".join(locations))
    keywords = uniprot_data.get("keywords", [])
    if keywords:
        lines.append("Keywords: " + ", ".join(keywords))
    go_terms = uniprot_data.get("go_terms", [])
    if go_terms:
        rendered = ", ".join(
            f"{t.get('aspect', '')}:{t.get('term', '')}" for t in go_terms
        )
        lines.append("GO terms: " + rendered)
    return "\n".join(lines)


async def extract_from_uniprot(gene: str, uniprot_data: dict[str, Any]) -> dict[str, Any]:
    """Format a UniProt record and extract annotations from it (no PMIDs)."""
    text = _format_uniprot_text(uniprot_data)
    return await extract_from_text(gene, text, source_pmids=[])


def _format_opentargets_text(ot_data: dict[str, Any]) -> str:
    """Render an OpenTargets record dict as a readable text block."""
    lines = [
        f"Ensembl ID: {ot_data.get('ensembl_id', '')}",
        f"Symbol: {ot_data.get('symbol', '')}",
        f"Approved name: {ot_data.get('approved_name', '')}",
    ]
    func_descriptions = ot_data.get("function_descriptions", [])
    if func_descriptions:
        lines.append("Function descriptions:")
        lines.extend(f"  - {d}" for d in func_descriptions)

    pathways = ot_data.get("pathways", [])
    if pathways:
        lines.append("Pathways:")
        lines.extend(
            f"  - {p.get('pathway', '')} ({p.get('pathwayId', '')})" for p in pathways
        )

    diseases = ot_data.get("disease_associations", [])
    if diseases:
        lines.append("Disease associations:")
        lines.extend(
            f"  - {d.get('name', '')} (score={d.get('score')})" for d in diseases
        )
    return "\n".join(lines)


async def extract_from_opentargets(gene: str, ot_data: dict[str, Any]) -> dict[str, Any]:
    """Format an OpenTargets record and extract annotations from it (no PMIDs)."""
    text = _format_opentargets_text(ot_data)
    return await extract_from_text(gene, text, source_pmids=[])
