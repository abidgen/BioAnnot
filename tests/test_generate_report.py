"""Tests for the self-contained HTML report generator (scripts/generate_report.py).

Writes a minimal run directory (final_annotations.json + prioritized_targets.tsv,
optionally plots/ and the cytoscape export) into tmp_path, renders the report, and
asserts the top-level tabs, per-gene content, priority table, network payload,
plots (base64 vs. placeholder), and pathway matrix.
"""

import base64
import json

import pytest

from scripts.generate_report import generate_report

# A 1×1 transparent PNG — enough to exercise base64 embedding of a real image.
_PNG_1x1 = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk"
    "+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=="
)


def _write_run(run_dir, annotations, tsv_rows, plots=None, cytoscape=None):
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "final_annotations.json").write_text(
        json.dumps(annotations), encoding="utf-8"
    )
    header = "gene\tcomposite\tsafety_flag\tconfidence\n"
    body = "".join(
        f"{g}\t{comp}\t{annotations[g]['safety_assessment'].get('safety_flag', False)}"
        f"\t{annotations[g]['confidence']}\n"
        for g, comp in tsv_rows
    )
    (run_dir / "prioritized_targets.tsv").write_text(header + body, encoding="utf-8")

    if plots:
        plots_dir = run_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        for name in plots:
            (plots_dir / name).write_bytes(_PNG_1x1)
    if cytoscape is not None:
        (run_dir / "target_network_cytoscape.json").write_text(
            json.dumps(cytoscape), encoding="utf-8"
        )


def _cytoscape_fixture():
    return {
        "elements": {
            "nodes": [
                {"data": {"id": "MYCN", "node_type": "target"}},
                {"data": {"id": "CD276", "node_type": "target"}},
                {"data": {"id": "MDM2", "node_type": "interactor"}},
            ],
            "edges": [
                {"data": {"source": "MYCN", "target": "CD276",
                          "type": "pathway_comembership"}},
                {"data": {"source": "MYCN", "target": "CD276",
                          "type": "direct_interaction"}},
                {"data": {"source": "MYCN", "target": "MDM2",
                          "type": "string_interaction"}},
            ],
        }
    }


@pytest.fixture
def annotations():
    return {
        "MYCN": {
            "gene_symbol": "MYCN",
            "functions": ["transcription activator"],
            "pathways": ["Signaling by WNT", "NON-CANONICAL: Ferroptosis"],
            "cellular_states": ["nucleus", "CellxGene: chromaffin cell"],
            "disease_associations": [
                {"disease": "neuroblastoma", "role": "oncogene",
                 "evidence_strength": "strong"},
            ],
            "interactors": ["USF1", "CDK2"],
            "druggability_notes": "fadraciclib shows response",
            "confidence": 0.84,
            "source_pmids": ["42233020", "42255978"],
            "source_count": 2,
            "merged_at": "2026-07-06T16:27:53Z",
            "safety_assessment": {
                "safety_flag": True,
                "tier1_flag": True,
                "tier2_flag": False,
                "tier1_high_tissues": {"Liver": 42.5},
                "safety_penalty": 0.60,
            },
            "cellxgene_expression": {
                "cell_type_count": 44, "tissue": "adrenal gland",
                "census_version": "2024-07-01",
                "top_cell_types": [{"cell_type": "chromaffin cell", "mean_expr": 0.0498}],
            },
        },
        "CD276": {
            "gene_symbol": "CD276",
            "functions": ["immune checkpoint"],
            "pathways": ["Signaling by WNT", "Immune System"],
            "cellular_states": ["cell surface"],
            "disease_associations": [
                {"disease": "neuroblastoma", "role": "therapeutic_target",
                 "evidence_strength": "moderate"},
            ],
            "interactors": [],
            "druggability_notes": "",
            "confidence": 0.70,
            "source_pmids": ["42270195"],
            "source_count": 1,
            "merged_at": "2026-07-06T16:28:10Z",
            "safety_assessment": {"safety_flag": False},
            "cellxgene_expression": {},
        },
    }


def _render(tmp_path, annotations, order=None, **kw):
    order = order or [("MYCN", 0.65), ("CD276", 0.47)]
    run_dir = tmp_path / "runs" / "20260706_122549"
    _write_run(run_dir, annotations, order, **kw)
    return generate_report(run_dir).read_text(encoding="utf-8")


def test_report_written_and_self_contained(tmp_path, annotations):
    html = _render(tmp_path, annotations, cytoscape=_cytoscape_fixture())
    assert "BioAnnot Annotation Report" in html
    assert "<link" not in html
    # The ONLY external reference is the D3 CDN script; data is inlined.
    assert html.count("<script src=") == 1
    assert "cdnjs.cloudflare.com/ajax/libs/d3" in html
    assert "window.__NETWORK__ =" in html


def test_top_level_tabs_render(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    for sec in ("genes", "priority", "network", "plots", "matrix"):
        assert f"data-sec='{sec}'" in html
        assert f"id='sec-{sec}'" in html
    for label in ("Priority Ranking", "Network Graph", "Plots", "Pathway Matrix"):
        assert label in html


def test_report_contains_gene_symbols_and_count(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    assert "MYCN" in html and "CD276" in html
    assert "Genes: <b>2</b>" in html


def test_report_has_pmid_and_genecards_links(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    assert "https://pubmed.ncbi.nlm.nih.gov/42233020/" in html
    assert "https://pubmed.ncbi.nlm.nih.gov/42270195/" in html
    assert "https://www.genecards.org/cgi-bin/carddisp.pl?gene=USF1" in html


def test_report_tier1_safety_banner(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    assert "VITAL ORGAN EXPRESSION" in html
    assert "Liver: 42.5 TPM" in html
    assert "safety-tier1" in html


def test_report_noncanonical_and_cellxgene_styling(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    assert "class='noncanon'" in html and "Ferroptosis" in html
    assert "class='cellxgene'" in html and "chromaffin cell" in html


def test_priority_ranking_table_sortable_with_tiers(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    assert "rank-table" in html
    assert "onclick='sortTable(this)'" in html
    # MYCN is tier-1 (red cell), CD276 clean (green cell).
    assert "tier-tier1" in html and "tier-clean" in html
    for col in ("Betweenness", "Degree", "Disease Score", "CellxGene Score", "Penalty"):
        assert col in html


def test_network_payload_inlined_with_edge_colors(tmp_path, annotations):
    html = _render(tmp_path, annotations, cytoscape=_cytoscape_fixture())
    # Data embedded (not fetched), so it survives a file:// open.
    payload = json.loads(
        html.split("window.__NETWORK__ = ", 1)[1].split(";</script>", 1)[0]
    )
    ids = {n["id"] for n in payload["nodes"]}
    assert {"MYCN", "CD276", "MDM2"} <= ids
    edge_colors = {l["color"] for l in payload["links"]}
    assert {"#028090", "#0D2B55", "#94A3B8"} <= edge_colors
    # Node color reflects safety tier (MYCN tier-1 red).
    mycn = next(n for n in payload["nodes"] if n["id"] == "MYCN")
    assert mycn["color"] == "#DC2626" and mycn["kind"] == "target"


def test_network_placeholder_when_cytoscape_absent(tmp_path, annotations):
    html = _render(tmp_path, annotations)  # no cytoscape file
    assert "window.__NETWORK__ = null" in html
    assert "export_cytoscape.py" in html


def test_plots_embedded_as_base64(tmp_path, annotations):
    html = _render(
        tmp_path, annotations,
        plots=["target_network.png", "score_breakdown.png",
               "pathway_heatmap.png", "cellxgene_expression.png"],
    )
    assert "data:image/png;base64," in html
    assert "Target Interaction Network" in html
    assert "Score Breakdown by Gene" in html
    # All four present → no placeholder text.
    assert "to generate this plot" not in html


def test_plot_placeholders_when_pngs_absent(tmp_path, annotations):
    html = _render(tmp_path, annotations)  # no plots/ dir
    # One placeholder per known plot (4).
    assert html.count("to generate this plot") == 4
    assert "visualize_network.py" in html


def test_pathway_matrix_excludes_noncanonical(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    matrix = html.split("id='sec-matrix'", 1)[1].split("</section>", 1)[0]
    assert "matrix-wrap" in matrix
    assert "Signaling by WNT" in matrix   # canonical, shared by both genes
    assert "Immune System" in matrix
    assert "Ferroptosis" not in matrix    # NON-CANONICAL excluded


def test_hash_routing_and_gene_selection(tmp_path, annotations):
    html = _render(tmp_path, annotations)
    assert "applyHash" in html and "hashchange" in html
    assert "selectGene('MYCN')" in html
    assert "data-gene='MYCN'" in html


def test_report_orders_genes_by_composite_desc(tmp_path, annotations):
    # Give CD276 the higher composite so order flips vs. dict order.
    html = _render(tmp_path, annotations, order=[("CD276", 0.90), ("MYCN", 0.20)])
    assert html.index("selectGene('CD276')") < html.index("selectGene('MYCN')")
    assert "class='tab active' data-gene='CD276'" in html
