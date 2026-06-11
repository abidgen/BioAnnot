"""Core LLM extraction module — Anthropic tool-use forced structured output.

Rules (per CLAUDE.md Step 5):
- Always use the extraction model and force the annotate_target tool.
- Never hallucinate PMIDs — only pass through IDs the fetcher returned.
- Truncate input text to ~3000 words before sending.
- Log token usage for every API call.
"""

from __future__ import annotations

import logging

import anthropic

from src.utils import validate_pmids

log = logging.getLogger("bio_annot.extractor")

EXTRACTION_MODEL = "claude-opus-4-8"

MAX_INPUT_WORDS = 3000
MAX_TOKENS = 4096

SYSTEM_PROMPT = (
    "You are a senior biomedical annotation scientist at NCI. Extract precise, "
    "evidence-based annotations only from the provided text. "
    "Use canonical gene symbols (HGNC). "
    "Use Reactome pathway names where possible. "
    "Set confidence < 0.5 if the text is tangentially related to the gene. "
    "Never invent interactors or pathways not mentioned in the text."
)

ANNOTATION_TOOL = [
    {
        "name": "annotate_target",
        "description": (
            "Extract structured biological annotations for a gene/protein from text"
        ),
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

_client: anthropic.Anthropic | None = None


def _get_client() -> anthropic.Anthropic:
    """Lazily construct the Anthropic client (resolves ANTHROPIC_API_KEY).

    Deferred so importing this module does not require the key to be set.
    """
    global _client
    if _client is None:
        _client = anthropic.Anthropic()
    return _client


def _truncate_words(text: str, max_words: int = MAX_INPUT_WORDS) -> str:
    """Truncate text to at most max_words whitespace-delimited words."""
    words = text.split()
    if len(words) <= max_words:
        return text
    log.info("Truncating input from %d to %d words", len(words), max_words)
    return " ".join(words[:max_words])


def extract_from_text(gene: str, text: str, source_pmids: list[str]) -> dict:
    """Extract structured annotations for a gene from free text.

    Calls the Anthropic API with the annotate_target tool forced, parses the
    tool_use block, attaches validated source PMIDs, and logs token usage.
    """
    truncated = _truncate_words(text)

    user_content = f"Gene: {gene}\n\n{truncated}"

    response = _get_client().messages.create(
        model=EXTRACTION_MODEL,
        max_tokens=MAX_TOKENS,
        system=SYSTEM_PROMPT,
        tools=ANNOTATION_TOOL,
        tool_choice={"type": "tool", "name": "annotate_target"},
        messages=[{"role": "user", "content": user_content}],
    )

    usage = response.usage
    log.info(
        "extract_from_text(%s): input_tokens=%s output_tokens=%s",
        gene,
        usage.input_tokens,
        usage.output_tokens,
    )

    tool_use = next(
        (block for block in response.content if block.type == "tool_use"), None
    )
    if tool_use is None:
        log.warning("No tool_use block returned for %s", gene)
        return {"gene_symbol": gene, "functions": [], "pathways": [], "confidence": 0.0}

    annotation = dict(tool_use.input)
    # Only pass through validated PMIDs — never invent or forward unvalidated IDs.
    annotation["source_pmids"] = validate_pmids(source_pmids)
    return annotation


def _format_uniprot_text(uniprot_data: dict) -> str:
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


def extract_from_uniprot(gene: str, uniprot_data: dict) -> dict:
    """Format a UniProt record and extract annotations from it (no PMIDs)."""
    text = _format_uniprot_text(uniprot_data)
    return extract_from_text(gene, text, source_pmids=[])


def _format_opentargets_text(ot_data: dict) -> str:
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


def extract_from_opentargets(gene: str, ot_data: dict) -> dict:
    """Format an OpenTargets record and extract annotations from it (no PMIDs)."""
    text = _format_opentargets_text(ot_data)
    return extract_from_text(gene, text, source_pmids=[])
