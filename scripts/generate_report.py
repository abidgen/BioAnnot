"""Generate a self-contained HTML annotation report for a pipeline run.

Reads ``final_annotations.json`` and ``prioritized_targets.tsv`` from the resolved
run directory (``outputs/latest`` by default; set ``RUN_DIR`` to target a specific
run) and writes ``bioannot_report.html`` alongside them. The report is a single
file with all CSS/JS inlined — no external assets, no server — so it opens by
double-clicking in any browser.

Run standalone::

    python scripts/generate_report.py

or import and call ``generate_report(run_dir)`` (the pipeline does this
automatically at the end of every run).
"""

from __future__ import annotations

import base64
import csv
import html
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

# Make the project root importable when run as a standalone script.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import load_disease_context, resolve_run_dir

# --- Color scheme (kept in one place; emitted as CSS custom properties) ---
COLORS = {
    "bg": "#F8FAFC",
    "navy": "#0D2B55",
    "teal": "#028090",
    "strong": "#02C39A",
    "moderate": "#F0A500",
    "weak": "#94A3B8",
    "tier1": "#DC2626",
    "tier2": "#D97706",
    "noncanon": "#F0A500",
}

# Disease-role → pill color.
ROLE_COLORS = {
    "oncogene": "#DC2626",
    "tumor_suppressor": "#028090",
    "biomarker": "#F0A500",
    "therapeutic_target": "#7C3AED",
    "unknown": "#94A3B8",
}

EVIDENCE_COLORS = {
    "strong": COLORS["strong"],
    "moderate": COLORS["moderate"],
    "weak": COLORS["weak"],
}

GENECARDS_URL = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}"
PUBMED_URL = "https://pubmed.ncbi.nlm.nih.gov/{pmid}/"


def _esc(value: Any) -> str:
    """HTML-escape any value's string form."""
    return html.escape("" if value is None else str(value))


def _slug(symbol: str) -> str:
    """DOM-safe id fragment for a gene symbol."""
    return "gene-" + re.sub(r"\W+", "_", symbol)


def _confidence_class(conf: float) -> str:
    """Grade a confidence value for badge coloring."""
    if conf >= 0.8:
        return "conf-high"
    if conf >= 0.6:
        return "conf-mid"
    return "conf-low"


def _read_annotations(run_dir: Path) -> dict[str, dict]:
    path = run_dir / "final_annotations.json"
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _read_cytoscape(run_dir: Path) -> dict | None:
    """Load the cytoscape network export, or None if absent/unreadable.

    Embedded inline into the report at build time (rather than fetched at
    runtime) so the network data survives a file:// double-click open.
    """
    path = run_dir / "target_network_cytoscape.json"
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def _read_priority_rows(
    run_dir: Path,
) -> tuple[list[str], dict[str, dict], dict[str, float]]:
    """Return (order by composite desc, {gene: full TSV row}, {gene: composite}).

    Missing/unreadable TSV yields empties; callers fall back to the annotation
    keys so the report still renders.
    """
    path = run_dir / "prioritized_targets.tsv"
    order: list[str] = []
    rows: dict[str, dict] = {}
    composite: dict[str, float] = {}
    if not path.exists():
        return order, rows, composite
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            order.append(gene)
            rows[gene] = row
            try:
                composite[gene] = float(row.get("composite", "") or "nan")
            except (ValueError, TypeError):
                composite[gene] = float("nan")
    return order, rows, composite


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in ("true", "1", "yes")


def _safety_tier(ann: dict, row: dict | None) -> tuple[str, str]:
    """Resolve (label, css-slug) safety tier from the annotation or TSV row.

    css-slug is one of 'tier1' / 'tier2' / 'clean'. Legacy records that are
    flagged but carry no tier detail are treated as tier 2 (amber), matching the
    per-gene safety banner.
    """
    sa = ann.get("safety_assessment") or {}
    row = row or {}
    tier1 = bool(sa.get("tier1_flag")) or _truthy(row.get("tier1_flag", ""))
    tier2 = bool(sa.get("tier2_flag")) or _truthy(row.get("tier2_flag", ""))
    flagged = bool(sa.get("safety_flag")) or _truthy(row.get("safety_flag", ""))
    if tier1:
        return "Tier 1", "tier1"
    if tier2:
        return "Tier 2", "tier2"
    if flagged:
        return "Tier 2", "tier2"
    return "Clean", "clean"


def _safety_penalty(ann: dict, row: dict | None, tier_css: str) -> float:
    """Resolve the composite penalty from the annotation/TSV, else the tier default."""
    sa = ann.get("safety_assessment") or {}
    if "safety_penalty" in sa:
        try:
            return float(sa["safety_penalty"])
        except (ValueError, TypeError):
            pass
    if row and row.get("safety_penalty"):
        try:
            return float(row["safety_penalty"])
        except (ValueError, TypeError):
            pass
    return {"tier1": 0.60, "tier2": 0.80, "clean": 1.00}[tier_css]


def _ordered_genes(
    annotations: dict[str, dict], tsv_order: list[str]
) -> list[str]:
    """Genes in TSV (composite-desc) order, with any TSV-absent genes appended."""
    seen = set()
    ordered = [g for g in tsv_order if g in annotations and not (g in seen or seen.add(g))]
    leftovers = sorted(
        (g for g in annotations if g not in seen),
        key=lambda g: annotations[g].get("confidence", 0.0),
        reverse=True,
    )
    return ordered + leftovers


def _run_timestamp(run_dir: Path, annotations: dict[str, dict]) -> str:
    """Human-readable run time from the run-dir name, else the latest merged_at."""
    name = run_dir.resolve().name
    try:
        return datetime.strptime(name, "%Y%m%d_%H%M%S").strftime("%Y-%m-%d %H:%M:%S")
    except ValueError:
        merged = [a.get("merged_at") for a in annotations.values() if a.get("merged_at")]
        return max(merged) if merged else "unknown"


# --- Section renderers (each returns an HTML fragment) ---

def _render_safety_banner(gene: dict) -> str:
    sa = gene.get("safety_assessment") or {}
    if not sa.get("safety_flag"):
        return ""

    if sa.get("tier1_flag"):
        cls, label = "safety-tier1", "VITAL ORGAN EXPRESSION"
        tissues = sa.get("tier1_high_tissues") or {}
    elif sa.get("tier2_flag"):
        cls, label = "safety-tier2", "SECONDARY TISSUE EXPRESSION"
        tissues = sa.get("tier2_high_tissues") or {}
    else:
        # Legacy pre-two-tier record: flagged but no tier detail.
        cls, label = "safety-tier2", "ELEVATED NORMAL-TISSUE EXPRESSION"
        tissues = {t: None for t in sa.get("high_expression_tissues") or []}

    items = []
    for tissue, tpm in tissues.items():
        if isinstance(tpm, (int, float)):
            items.append(f"<span class='tissue'>{_esc(tissue)}: {tpm:.1f} TPM</span>")
        else:
            items.append(f"<span class='tissue'>{_esc(tissue)}</span>")
    body = " ".join(items) if items else "<span class='tissue'>(tissues unspecified)</span>"
    return (
        f"<div class='safety-banner {cls}'>"
        f"<span class='safety-label'>⚠ {label}</span>{body}</div>"
    )


def _render_list(items: list, item_class_fn=None) -> str:
    if not items:
        return "<p class='empty'>None reported.</p>"
    lis = []
    for item in items:
        cls = item_class_fn(item) if item_class_fn else ""
        text = _esc(item)
        if cls == "noncanon":
            text = "⚠ " + text
        cls_attr = f" class='{cls}'" if cls else ""
        lis.append(f"<li{cls_attr}>{text}</li>")
    return "<ul>" + "".join(lis) + "</ul>"


def _cellular_state_class(state: str) -> str:
    return "cellxgene" if str(state).startswith("CellxGene:") else ""


def _pathway_class(pathway: str) -> str:
    return "noncanon" if str(pathway).startswith("NON-CANONICAL:") else ""


def _render_disease_table(assocs: list[dict]) -> str:
    if not assocs:
        return "<p class='empty'>No disease associations.</p>"
    rows = []
    for a in assocs:
        role = str(a.get("role", "unknown"))
        color = ROLE_COLORS.get(role, ROLE_COLORS["unknown"])
        strength = str(a.get("evidence_strength", "weak")).lower()
        pill = (
            f"<span class='pill' style='background:{color}'>"
            f"{_esc(role.replace('_', ' '))}</span>"
        )
        rows.append(
            f"<tr class='ev-{_esc(strength)}'>"
            f"<td>{_esc(a.get('disease', ''))}</td>"
            f"<td>{pill}</td>"
            f"<td class='ev-label'>{_esc(strength.upper())}</td></tr>"
        )
    return (
        "<table class='disease-table'><thead><tr>"
        "<th>Disease</th><th>Role</th><th>Evidence</th></tr></thead>"
        "<tbody>" + "".join(rows) + "</tbody></table>"
    )


def _render_interactors(interactors: list) -> str:
    if not interactors:
        return "<p class='empty'>None reported.</p>"
    tags = []
    for sym in interactors:
        url = GENECARDS_URL.format(symbol=_esc(sym))
        tags.append(
            f"<a class='tag' href='{url}' target='_blank' rel='noopener'>{_esc(sym)}</a>"
        )
    return "<div class='tags'>" + "".join(tags) + "</div>"


def _render_pmids(pmids: list) -> str:
    if not pmids:
        return "<p class='empty'>No source PMIDs.</p>"
    links = []
    for pmid in pmids:
        url = PUBMED_URL.format(pmid=_esc(pmid))
        links.append(
            f"<a class='pmid' href='{url}' target='_blank' rel='noopener'>{_esc(pmid)}</a>"
        )
    return "<div class='pmids'>" + "".join(links) + "</div>"


def _render_cellxgene(cx: dict) -> str:
    if not cx:
        return ""
    top = cx.get("top_cell_types") or []
    if not top:
        return ""
    rows = "".join(
        f"<tr><td>{_esc(c.get('cell_type', ''))}</td>"
        f"<td>{float(c.get('mean_expr', 0.0)):.4f}</td></tr>"
        for c in top
    )
    caption = _esc(
        f"tissue: {cx.get('tissue', '?')} · {cx.get('cell_type_count', '?')} cell types "
        f"· census {cx.get('census_version', '?')}"
    )
    return (
        "<div class='section'><h3>CellxGene Expression</h3>"
        f"<p class='caption'>{caption}</p>"
        "<table class='cx-table'><thead><tr><th>Cell Type</th>"
        "<th>Mean Expr</th></tr></thead><tbody>" + rows + "</tbody></table></div>"
    )


def _render_druggability(notes: str) -> str:
    if not notes:
        return ""
    return (
        "<div class='section'><h3>Druggability Notes</h3>"
        f"<p class='druggability'>{_esc(notes)}</p></div>"
    )


def _render_gene_panel(symbol: str, gene: dict, composite: float, active: bool) -> str:
    conf = float(gene.get("confidence", 0.0) or 0.0)
    comp_txt = "—" if composite != composite else f"{composite:.3f}"  # NaN check
    header = (
        f"<div class='gene-head'>"
        f"<h2 class='gene-symbol'>{_esc(symbol)}</h2>"
        f"<span class='badge {_confidence_class(conf)}'>confidence {conf:.2f}</span>"
        f"<span class='meta'>sources: {_esc(gene.get('source_count', '?'))}</span>"
        f"<span class='meta'>composite: {comp_txt}</span>"
        f"<span class='meta'>merged: {_esc(gene.get('merged_at', 'unknown'))}</span>"
        f"</div>"
    )

    left = (
        "<div class='col'>"
        "<div class='section'><h3>Functions</h3>"
        f"{_render_list(gene.get('functions', []))}</div>"
        "<div class='section'><h3>Cellular States</h3>"
        f"{_render_list(gene.get('cellular_states', []), _cellular_state_class)}</div>"
        "</div>"
    )
    right = (
        "<div class='col'>"
        "<div class='section'><h3>Pathways</h3>"
        f"{_render_list(gene.get('pathways', []), _pathway_class)}</div>"
        f"{_render_druggability(gene.get('druggability_notes', ''))}"
        "</div>"
    )

    body = (
        header
        + _render_safety_banner(gene)
        + f"<div class='columns'>{left}{right}</div>"
        + "<div class='section'><h3>Disease Associations</h3>"
        + _render_disease_table(gene.get("disease_associations", []))
        + "</div>"
        + "<div class='section'><h3>Interactors</h3>"
        + _render_interactors(gene.get("interactors", []))
        + "</div>"
        + "<div class='section'><h3>Sources</h3>"
        + _render_pmids(gene.get("source_pmids", []))
        + "</div>"
        + _render_cellxgene(gene.get("cellxgene_expression") or {})
    )

    active_cls = " active" if active else ""
    return (
        f"<section id='{_slug(symbol)}' class='gene-panel{active_cls}' "
        f"data-gene='{_esc(symbol)}'>{body}</section>"
    )


# Safety-tier → color (shared by ranking cells and network nodes).
TIER_HEX = {
    "tier1": COLORS["tier1"],
    "tier2": COLORS["tier2"],
    "clean": COLORS["strong"],
}
EDGE_HEX = {
    "pathway_comembership": COLORS["teal"],
    "direct_interaction": COLORS["navy"],
    "string_interaction": COLORS["weak"],
}

# Known plots (filename → card title), in display order.
PLOT_SPECS = [
    ("target_network.png", "Target Interaction Network"),
    ("score_breakdown.png", "Score Breakdown by Gene"),
    ("pathway_heatmap.png", "Pathway Co-membership Heatmap"),
    ("cellxgene_expression.png", "CellxGene Single-Cell Expression"),
]


def _fmt(value: Any, digits: int = 3) -> str:
    """Format a numeric cell to N decimals, or '—' if empty/non-numeric."""
    if value is None or value == "":
        return "—"
    try:
        return f"{float(value):.{digits}f}"
    except (ValueError, TypeError):
        return _esc(value)


def _render_priority_table(
    genes: list[str], annotations: dict, rows: dict, composite: dict
) -> str:
    """Sortable ranking table (one row per gene, composite-desc by default)."""
    cols = [
        ("Rank", "num"), ("Gene", "str"), ("Composite", "num"),
        ("Betweenness", "num"), ("Degree", "num"), ("Disease Score", "num"),
        ("CellxGene Score", "num"), ("Confidence", "num"),
        ("Druggability", "num"), ("Safety Tier", "str"), ("Penalty", "num"),
    ]
    head = "".join(
        f"<th data-type='{t}' onclick='sortTable(this)'>{_esc(label)}"
        f"<span class='arrow'></span></th>"
        for label, t in cols
    )
    body_rows = []
    for i, g in enumerate(genes, start=1):
        row = rows.get(g, {})
        ann = annotations.get(g, {})
        label, css = _safety_tier(ann, row)
        penalty = _safety_penalty(ann, row, css)
        comp = composite.get(g, float("nan"))
        comp_txt = "—" if comp != comp else f"{comp:.3f}"
        cells = [
            f"<td data-v='{i}'>{i}</td>",
            f"<td data-v='{_esc(g)}'>{_esc(g)}</td>",
            f"<td data-v='{comp if comp == comp else -1}'>{comp_txt}</td>",
            f"<td data-v='{_esc(row.get('betweenness',''))}'>{_fmt(row.get('betweenness'))}</td>",
            f"<td data-v='{_esc(row.get('degree',''))}'>{_fmt(row.get('degree'))}</td>",
            f"<td data-v='{_esc(row.get('disease_score',''))}'>{_fmt(row.get('disease_score'),2)}</td>",
            f"<td data-v='{_esc(row.get('cellxgene_score',''))}'>{_fmt(row.get('cellxgene_score'),2)}</td>",
            f"<td data-v='{_esc(row.get('confidence',''))}'>{_fmt(row.get('confidence'),2)}</td>",
            f"<td data-v='{_esc(row.get('druggability_bonus',''))}'>{_fmt(row.get('druggability_bonus'),2)}</td>",
            f"<td data-v='{_esc(label)}' class='tier-cell tier-{css}'>{_esc(label)}</td>",
            f"<td data-v='{penalty:.2f}'>{penalty:.2f}</td>",
        ]
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        "<table class='rank-table'><thead><tr>" + head + "</tr></thead>"
        "<tbody>" + "".join(body_rows) + "</tbody></table>"
    )


def _build_network_payload(
    cyto: dict | None, annotations: dict, rows: dict, composite: dict
) -> dict | None:
    """Build a compact {nodes, links} payload for D3 from the cytoscape export."""
    if not cyto:
        return None
    elements = cyto.get("elements", {})
    nodes = []
    for n in elements.get("nodes", []):
        d = n.get("data", {})
        gid = d.get("id")
        if gid is None:
            continue
        if d.get("node_type") == "target":
            _, css = _safety_tier(annotations.get(gid, {}), rows.get(gid))
            comp = composite.get(gid)
            comp = None if (comp is None or comp != comp) else comp
            nodes.append({
                "id": gid, "kind": "target", "composite": comp,
                "tier": {"tier1": "Tier 1", "tier2": "Tier 2", "clean": "Clean"}[css],
                "color": TIER_HEX[css],
            })
        else:
            nodes.append({
                "id": gid, "kind": "interactor", "composite": None,
                "tier": "interactor", "color": COLORS["weak"],
            })
    links = [
        {
            "source": e["data"].get("source"),
            "target": e["data"].get("target"),
            "color": EDGE_HEX.get(e["data"].get("type"), "#cccccc"),
            "type": e["data"].get("type"),
        }
        for e in elements.get("edges", [])
        if e.get("data")
    ]
    return {"nodes": nodes, "links": links}


def _render_network_section(payload: dict | None) -> str:
    if not payload:
        return (
            "<div class='placeholder'>Run <code>export_cytoscape.py</code> to "
            "generate the network data (target_network_cytoscape.json).</div>"
        )
    legend = (
        "<div class='legend'>"
        f"<span><i style='background:{COLORS['tier1']}'></i>Tier 1</span>"
        f"<span><i style='background:{COLORS['tier2']}'></i>Tier 2</span>"
        f"<span><i style='background:{COLORS['strong']}'></i>Clean</span>"
        f"<span><i style='background:{COLORS['weak']}'></i>Interactor</span>"
        f"<span><i class='edge' style='background:{COLORS['teal']}'></i>pathway</span>"
        f"<span><i class='edge' style='background:{COLORS['navy']}'></i>interaction</span>"
        f"<span><i class='edge' style='background:{COLORS['weak']}'></i>STRING</span>"
        "</div>"
    )
    return (
        "<p class='caption'>Node size ∝ composite score · color = safety tier · "
        "click a target to open its annotation · scroll to zoom, drag to pan.</p>"
        + legend
        + "<div class='net-controls'>"
        "<button id='net-reset' onclick='resetNetworkView()'>Reset View</button>"
        "</div>"
        "<div id='net-wrap'>"
        "<svg id='net-svg' width='100%' height='100%' "
        "viewBox='0 0 960 620' preserveAspectRatio='xMidYMid meet'></svg>"
        "<div id='net-tip'></div></div>"
    )


def _render_plots_section(run_dir: Path) -> str:
    cards = []
    plots_dir = run_dir / "plots"
    for filename, title in PLOT_SPECS:
        path = plots_dir / filename
        if path.exists():
            b64 = base64.b64encode(path.read_bytes()).decode("ascii")
            inner = f"<img src='data:image/png;base64,{b64}' alt='{_esc(title)}'>"
        else:
            inner = (
                "<div class='placeholder'>Run <code>visualize_network.py</code> "
                "to generate this plot</div>"
            )
        cards.append(
            f"<div class='plot-card'><h3>{_esc(title)}</h3>{inner}</div>"
        )
    return "<div class='plot-stack'>" + "".join(cards) + "</div>"


def _render_pathway_matrix(genes: list[str], annotations: dict) -> str:
    """genes × canonical-pathway presence matrix (rows composite-desc, cols by freq)."""
    gene_paths = {
        g: {
            p for p in (annotations.get(g, {}).get("pathways") or [])
            if not str(p).startswith("NON-CANONICAL:")
        }
        for g in genes
    }
    freq: Counter = Counter()
    for paths in gene_paths.values():
        freq.update(paths)
    if not freq:
        return "<p class='empty'>No canonical pathways to display.</p>"
    cols = sorted(freq, key=lambda p: (-freq[p], p))

    header = "<th class='corner'></th>" + "".join(
        f"<th class='colh'>{_esc(p)}</th>" for p in cols
    )
    body_rows = []
    for g in genes:
        cells = [f"<th class='rowh'>{_esc(g)}</th>"]
        for p in cols:
            if p in gene_paths[g]:
                cells.append(
                    f"<td class='cell on' title='{_esc(g)} — {_esc(p)}'></td>"
                )
            else:
                cells.append(f"<td class='cell' title='{_esc(g)} — {_esc(p)}'></td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    legend = (
        "<div class='matrix-legend'>"
        "<span><i class='sq on'></i>Pathway present</span>"
        "<span><i class='sq'></i>Pathway absent</span>"
        "</div>"
    )
    return (
        "<div class='matrix-wrap'><table class='matrix'><thead><tr>"
        + header
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table></div>"
        + legend
    )


_ROOT_VARS = ":root {" + "".join(f"--{k}:{v};" for k, v in COLORS.items()) + "}"

CSS = _ROOT_VARS + """
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--bg); color: #1E293B;
  font-family: system-ui, -apple-system, sans-serif; line-height: 1.5;
}
header.app {
  background: var(--navy); color: #fff; padding: 20px 28px;
}
header.app h1 { margin: 0 0 6px; font-size: 24px; }
header.app .sub { opacity: .85; font-size: 14px; }
header.app .sub b { color: #7FE0D6; }
nav.tabs {
  position: sticky; top: 0; z-index: 10; background: var(--navy);
  padding: 0 12px; display: flex; gap: 4px; overflow-x: auto;
  border-top: 1px solid rgba(255,255,255,.12);
}
nav.tabs .tab {
  background: transparent; color: #cbd5e1; border: 0; cursor: pointer;
  padding: 12px 16px; font-size: 14px; font-weight: 600; white-space: nowrap;
  border-bottom: 3px solid transparent;
}
nav.tabs .tab:hover { color: #fff; }
nav.tabs .tab.active { color: #fff; border-bottom-color: var(--teal); }
main { max-width: 1100px; margin: 0 auto; padding: 24px 20px 60px; }
.gene-panel { display: none; }
.gene-panel.active { display: block; }
.gene-head {
  display: flex; align-items: center; flex-wrap: wrap; gap: 12px;
  padding-bottom: 14px; border-bottom: 2px solid #E2E8F0; margin-bottom: 16px;
}
.gene-symbol { margin: 0; font-size: 34px; color: var(--navy); letter-spacing: .5px; }
.badge {
  color: #fff; padding: 4px 10px; border-radius: 999px; font-size: 12px;
  font-weight: 700;
}
.conf-high { background: var(--strong); }
.conf-mid { background: var(--moderate); }
.conf-low { background: var(--weak); }
.meta { font-size: 13px; color: #475569; }
.columns { display: flex; gap: 24px; flex-wrap: wrap; }
.col { flex: 1; min-width: 300px; }
.section { margin-bottom: 20px; }
.section h3 {
  font-size: 15px; text-transform: uppercase; letter-spacing: .6px;
  color: var(--teal); margin: 0 0 8px; border-bottom: 1px solid #E2E8F0;
  padding-bottom: 4px;
}
ul { margin: 0; padding-left: 20px; }
li { margin: 3px 0; }
li.cellxgene { color: var(--teal); font-style: italic; }
li.noncanon { color: var(--noncanon); font-weight: 600; }
.empty { color: #94A3B8; font-style: italic; margin: 0; }
.druggability {
  background: #fff; border: 1px solid #E2E8F0; border-radius: 8px;
  padding: 12px 14px; font-style: italic; margin: 0;
}
.safety-banner {
  display: flex; align-items: center; flex-wrap: wrap; gap: 10px;
  color: #fff; padding: 12px 16px; border-radius: 8px; margin-bottom: 16px;
  font-size: 14px;
}
.safety-banner.safety-tier1 { background: var(--tier1); }
.safety-banner.safety-tier2 { background: var(--tier2); }
.safety-label { font-weight: 800; letter-spacing: .5px; }
.safety-banner .tissue {
  background: rgba(255,255,255,.22); padding: 2px 8px; border-radius: 4px;
}
table { border-collapse: collapse; width: 100%; background: #fff; }
.disease-table th, .cx-table th {
  text-align: left; background: #F1F5F9; color: #334155; font-size: 12px;
  text-transform: uppercase; letter-spacing: .4px; padding: 8px 10px;
}
.disease-table td, .cx-table td {
  padding: 8px 10px; border-top: 1px solid #E2E8F0; font-size: 14px;
}
.disease-table tr.ev-strong td:first-child { border-left: 4px solid var(--strong); }
.disease-table tr.ev-moderate td:first-child { border-left: 4px solid var(--moderate); }
.disease-table tr.ev-weak td:first-child { border-left: 4px solid var(--weak); }
.ev-label { font-weight: 700; font-size: 12px; }
.pill {
  color: #fff; padding: 3px 10px; border-radius: 999px; font-size: 12px;
  font-weight: 700; white-space: nowrap;
}
.tags, .pmids { display: flex; flex-wrap: wrap; gap: 8px; }
.tag {
  background: #fff; border: 1px solid var(--teal); color: var(--teal);
  padding: 4px 12px; border-radius: 999px; font-size: 13px; font-weight: 600;
  text-decoration: none;
}
.tag:hover { background: var(--teal); color: #fff; }
.pmid {
  background: #EEF2FF; color: var(--navy); padding: 4px 10px; border-radius: 6px;
  font-size: 13px; text-decoration: none; font-variant-numeric: tabular-nums;
}
.pmid:hover { background: var(--navy); color: #fff; }
.caption { color: #64748B; font-size: 12px; margin: 0 0 8px; }
.cx-table { max-width: 480px; }

/* Top-level section navigation */
nav.toptabs {
  background: var(--navy); display: flex; gap: 6px; padding: 8px 16px 0;
  overflow-x: auto;
}
nav.toptabs .toptab {
  background: #E2E8F0; color: #475569; border: 0; cursor: pointer;
  padding: 10px 18px; font-size: 14px; font-weight: 700; white-space: nowrap;
  border-radius: 8px 8px 0 0;
}
nav.toptabs .toptab.active { background: var(--navy); color: #fff; }
.topsec { display: none; }
.topsec.active { display: block; }

/* Priority ranking table */
.rank-table { font-size: 13px; }
.rank-table th {
  background: var(--navy); color: #fff; text-align: right; padding: 8px 10px;
  cursor: pointer; user-select: none; white-space: nowrap; position: sticky; top: 0;
}
.rank-table th:nth-child(2), .rank-table th:nth-child(10) { text-align: left; }
.rank-table td { padding: 6px 10px; border-top: 1px solid #E2E8F0; text-align: right;
  font-variant-numeric: tabular-nums; }
.rank-table td:nth-child(2) { text-align: left; font-weight: 700; color: var(--navy); }
.rank-table tbody tr:hover { background: #EEF2FF; }
.rank-table .arrow { font-size: 10px; margin-left: 4px; opacity: .7; }
.tier-cell { color: #fff; text-align: center; font-weight: 700; }
.tier-tier1 { background: var(--tier1); }
.tier-tier2 { background: var(--tier2); }
.tier-clean { background: var(--strong); }

/* Network graph */
.placeholder {
  background: #fff; border: 2px dashed #CBD5E1; border-radius: 10px;
  padding: 40px; text-align: center; color: #64748B; font-size: 15px;
}
.placeholder code { background: #F1F5F9; padding: 2px 6px; border-radius: 4px; }
.legend { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 10px; font-size: 13px; }
.legend span { display: flex; align-items: center; gap: 6px; }
.legend i { width: 14px; height: 14px; border-radius: 50%; display: inline-block; }
.legend i.edge { width: 18px; height: 3px; border-radius: 0; }
.net-controls { margin-bottom: 8px; }
#net-reset {
  background: var(--navy); color: #fff; border: 0; cursor: pointer;
  padding: 6px 14px; border-radius: 6px; font-size: 13px; font-weight: 600;
}
#net-reset:hover { background: var(--teal); }
#net-wrap { position: relative; background: #fff; border: 1px solid #E2E8F0;
  border-radius: 10px; overflow: hidden;
  width: 100%; height: 70vh; min-height: 500px; }
#net-svg { width: 100%; height: 100%; display: block; cursor: grab; }
#net-svg:active { cursor: grabbing; }
#net-tip {
  position: absolute; pointer-events: none; background: var(--navy); color: #fff;
  padding: 6px 10px; border-radius: 6px; font-size: 12px; opacity: 0;
  transition: opacity .1s; white-space: nowrap; z-index: 5;
}

/* Plots */
.plot-stack { display: flex; flex-direction: column; }
.plot-card { background: #fff; border: 1px solid #E2E8F0; border-radius: 10px;
  padding: 14px; width: 100%; margin-bottom: 32px; }
.plot-card h3 { margin: 0 0 10px; color: var(--teal); font-size: 15px;
  text-transform: uppercase; letter-spacing: .5px; }
.plot-card img { width: 100%; max-width: 100%; height: auto; display: block;
  border-radius: 8px; }

/* Pathway matrix */
.matrix-wrap { overflow-x: auto; overflow-y: auto; max-height: 78vh;
  border: 1px solid #E2E8F0; border-radius: 10px; background: #fff; }
.matrix { border-collapse: collapse; }
.matrix th.corner { position: sticky; left: 0; top: 0; z-index: 3;
  background: var(--bg); }
.matrix th.rowh { position: sticky; left: 0; z-index: 1; background: var(--bg);
  text-align: left; font-weight: 700; padding-right: 12px; padding-left: 10px;
  white-space: nowrap; font-size: 13px; color: var(--navy);
  border-top: 1px solid #E2E8F0; }
.matrix th.colh {
  writing-mode: vertical-rl; transform: rotate(180deg); white-space: nowrap;
  height: 180px; max-height: 180px; font-size: 11px; font-weight: 600;
  vertical-align: bottom; padding: 6px 0; color: #334155;
  position: sticky; top: 0; background: #fff; z-index: 2;
  min-width: 28px;
}
.matrix td.cell { min-width: 28px; width: 28px; height: 28px;
  border: 1px solid #EEF2F6; }
.matrix td.cell.on { background: var(--teal); }
.matrix-legend { display: flex; gap: 18px; margin-top: 12px; font-size: 13px; }
.matrix-legend span { display: flex; align-items: center; gap: 6px; }
.matrix-legend i.sq { width: 16px; height: 16px; display: inline-block;
  border: 1px solid #CBD5E1; }
.matrix-legend i.sq.on { background: var(--teal); border-color: var(--teal); }
"""


JS = r"""
var SECTIONS = ['genes', 'priority', 'network', 'plots', 'matrix'];
var CURRENT_GENE = null;

function activateSection(sec) {
  document.querySelectorAll('.topsec').forEach(function (s) { s.classList.remove('active'); });
  document.querySelectorAll('nav.toptabs .toptab').forEach(function (t) { t.classList.remove('active'); });
  var el = document.getElementById('sec-' + sec);
  if (el) el.classList.add('active');
  var tab = document.querySelector(".toptab[data-sec='" + sec + "']");
  if (tab) tab.classList.add('active');
  if (sec === 'network') renderNetwork();
}

function activateGene(sym) {
  CURRENT_GENE = sym;
  document.querySelectorAll('.gene-panel').forEach(function (p) { p.classList.remove('active'); });
  document.querySelectorAll('nav.tabs .tab').forEach(function (t) { t.classList.remove('active'); });
  var panel = document.querySelector(".gene-panel[data-gene='" + sym + "']");
  if (panel) panel.classList.add('active');
  var tab = document.querySelector(".tab[data-gene='" + sym + "']");
  if (tab) tab.classList.add('active');
}

function selectGene(sym) { location.hash = 'genes/' + sym; }

function goSec(sec) {
  location.hash = (sec === 'genes') ? ('genes/' + (CURRENT_GENE || FIRST_GENE)) : sec;
}

function applyHash() {
  var raw = (location.hash || '').replace(/^#\/?/, '');
  var parts = raw.split('/');
  var sec = parts[0];
  if (SECTIONS.indexOf(sec) < 0) sec = 'genes';
  activateSection(sec);
  if (sec === 'genes') activateGene(parts[1] || FIRST_GENE);
}

function sortTable(th) {
  var table = th.closest('table');
  var idx = Array.prototype.indexOf.call(th.parentNode.children, th);
  var type = th.getAttribute('data-type') || 'str';
  var asc = th.getAttribute('data-asc') === 'true' ? false : true;
  Array.prototype.forEach.call(th.parentNode.children, function (h) {
    h.removeAttribute('data-asc');
    var a = h.querySelector('.arrow'); if (a) a.textContent = '';
  });
  th.setAttribute('data-asc', asc);
  var arrow = th.querySelector('.arrow'); if (arrow) arrow.textContent = asc ? ' ▲' : ' ▼';
  var tbody = table.querySelector('tbody');
  var rows = Array.prototype.slice.call(tbody.querySelectorAll('tr'));
  rows.sort(function (a, b) {
    var av = a.children[idx].getAttribute('data-v');
    var bv = b.children[idx].getAttribute('data-v');
    if (type === 'num') {
      av = parseFloat(av); bv = parseFloat(bv);
      if (isNaN(av)) av = -Infinity; if (isNaN(bv)) bv = -Infinity;
      return asc ? av - bv : bv - av;
    }
    av = (av || '').toLowerCase(); bv = (bv || '').toLowerCase();
    if (av < bv) return asc ? -1 : 1;
    if (av > bv) return asc ? 1 : -1;
    return 0;
  });
  rows.forEach(function (r) { tbody.appendChild(r); });
}

function resetNetworkView() {
  if (window.__net_svg && window.__net_zoom) {
    window.__net_svg.transition().duration(400)
      .call(window.__net_zoom.transform, d3.zoomIdentity);
  }
}

function renderNetwork() {
  if (window.__netDone) return;
  var data = window.__NETWORK__;
  if (!data || typeof d3 === 'undefined') return;
  window.__netDone = true;
  var W = 960, H = 620;
  var svg = d3.select('#net-svg');
  var g = svg.append('g');
  // Zoom/pan: scroll to zoom (0.3×–4×), drag background to pan. Stored on window
  // so the "Reset View" button can restore the identity transform.
  var zoom = d3.zoom().scaleExtent([0.3, 4])
    .on('zoom', function (e) { g.attr('transform', e.transform); });
  svg.call(zoom);
  window.__net_svg = svg;
  window.__net_zoom = zoom;
  function radius(d) { return d.composite != null ? 7 + d.composite * 44 : 5; }
  var link = g.append('g').attr('stroke-opacity', 0.55).selectAll('line')
    .data(data.links).join('line')
    .attr('stroke', function (d) { return d.color; }).attr('stroke-width', 1.3);
  var node = g.append('g').selectAll('circle').data(data.nodes).join('circle')
    .attr('r', radius).attr('fill', function (d) { return d.color; })
    .attr('stroke', '#fff').attr('stroke-width', 1.2)
    .style('cursor', function (d) { return d.kind === 'target' ? 'pointer' : 'default'; })
    .on('click', function (e, d) { if (d.kind === 'target') location.hash = 'genes/' + d.id; })
    .on('mouseover', function (e, d) { showTip(d); })
    .on('mousemove', moveTip)
    .on('mouseout', function () { document.getElementById('net-tip').style.opacity = 0; })
    .call(d3.drag()
      .on('start', function (e, d) { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag', function (e, d) { d.fx = e.x; d.fy = e.y; })
      .on('end', function (e, d) { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));
  var label = g.append('g').selectAll('text')
    .data(data.nodes.filter(function (d) { return d.kind === 'target'; })).join('text')
    .text(function (d) { return d.id; }).attr('font-size', 10).attr('fill', '#0D2B55')
    .attr('dx', 7).attr('dy', 3).style('pointer-events', 'none');
  var sim = d3.forceSimulation(data.nodes)
    .force('link', d3.forceLink(data.links).id(function (d) { return d.id; }).distance(72))
    .force('charge', d3.forceManyBody().strength(-170))
    .force('center', d3.forceCenter(W / 2, H / 2))
    .force('collide', d3.forceCollide().radius(function (d) { return radius(d) + 3; }))
    .on('tick', function () {
      link.attr('x1', function (d) { return d.source.x; }).attr('y1', function (d) { return d.source.y; })
        .attr('x2', function (d) { return d.target.x; }).attr('y2', function (d) { return d.target.y; });
      node.attr('cx', function (d) { return d.x; }).attr('cy', function (d) { return d.y; });
      label.attr('x', function (d) { return d.x; }).attr('y', function (d) { return d.y; });
    });
  var tipEl = document.getElementById('net-tip');
  function showTip(d) {
    tipEl.innerHTML = '<b>' + d.id + '</b><br>composite: ' +
      (d.composite != null ? d.composite.toFixed(3) : '—') + '<br>safety: ' + d.tier;
    tipEl.style.opacity = 1;
  }
  function moveTip(e) {
    var wrap = document.getElementById('net-wrap').getBoundingClientRect();
    tipEl.style.left = (e.clientX - wrap.left + 12) + 'px';
    tipEl.style.top = (e.clientY - wrap.top + 12) + 'px';
  }
}

window.addEventListener('hashchange', applyHash);
window.addEventListener('DOMContentLoaded', function () {
  if (!location.hash) { location.hash = 'genes/' + FIRST_GENE; } else { applyHash(); }
});
"""


D3_CDN = "https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"

TOP_TABS = [
    ("genes", "Genes"),
    ("priority", "Priority Ranking"),
    ("network", "Network Graph"),
    ("plots", "Plots"),
    ("matrix", "Pathway Matrix"),
]


def render_html(
    annotations: dict[str, dict],
    genes: list[str],
    rows: dict[str, dict],
    composite: dict[str, float],
    network_payload: dict | None,
    run_dir: Path,
    disease_context: str,
    run_ts: str,
) -> str:
    # --- Genes section (per-gene sub-tabs + panels) ---
    gene_tabs = "".join(
        f"<button class='tab{' active' if i == 0 else ''}' data-gene='{_esc(g)}' "
        f"onclick=\"selectGene('{_esc(g)}')\">{_esc(g)}</button>"
        for i, g in enumerate(genes)
    )
    gene_panels = "".join(
        _render_gene_panel(
            g, annotations[g], composite.get(g, float("nan")), active=(i == 0)
        )
        for i, g in enumerate(genes)
    )
    genes_section = (
        "<section class='topsec active' id='sec-genes'>"
        f"<nav class='tabs'>{gene_tabs}</nav>{gene_panels}</section>"
    )

    # --- Other sections ---
    priority_section = (
        "<section class='topsec' id='sec-priority'>"
        "<div class='section'><h3>Priority Ranking</h3>"
        + _render_priority_table(genes, annotations, rows, composite)
        + "</div></section>"
    )
    network_section = (
        "<section class='topsec' id='sec-network'>"
        "<div class='section'><h3>Network Graph</h3>"
        + _render_network_section(network_payload)
        + "</div></section>"
    )
    plots_section = (
        "<section class='topsec' id='sec-plots'>"
        "<div class='section'><h3>Plots</h3>"
        + _render_plots_section(run_dir)
        + "</div></section>"
    )
    matrix_section = (
        "<section class='topsec' id='sec-matrix'>"
        "<div class='section'><h3>Pathway Matrix</h3>"
        + _render_pathway_matrix(genes, annotations)
        + "</div></section>"
    )

    top_nav = "<nav class='toptabs'>" + "".join(
        f"<button class='toptab{' active' if sec == 'genes' else ''}' "
        f"data-sec='{sec}' onclick=\"goSec('{sec}')\">{_esc(label)}</button>"
        for sec, label in TOP_TABS
    ) + "</nav>"

    header = (
        "<header class='app'><h1>BioAnnot Annotation Report</h1>"
        f"<div class='sub'>Disease context: <b>{_esc(disease_context)}</b> &nbsp;·&nbsp; "
        f"Run: <b>{_esc(run_ts)}</b> &nbsp;·&nbsp; "
        f"Genes: <b>{len(genes)}</b></div></header>"
    )

    first_gene = genes[0] if genes else ""
    # Escape "<" so an embedded value can never break out of the <script> tag
    # (json.dumps does not escape forward slashes / "</script>").
    payload_js = json.dumps(network_payload).replace("<", "\\u003c")
    first_gene_js = json.dumps(first_gene).replace("<", "\\u003c")
    data_script = (
        f"<script>var FIRST_GENE = {first_gene_js}; "
        f"window.__NETWORK__ = {payload_js};</script>"
    )

    return (
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        "<title>BioAnnot Annotation Report</title>"
        f"<style>{CSS}</style></head><body>"
        f"{header}{top_nav}<main>"
        f"{genes_section}{priority_section}{network_section}{plots_section}{matrix_section}"
        "</main>"
        f"{data_script}"
        f"<script src='{D3_CDN}'></script>"
        f"<script>{JS}</script></body></html>"
    )


def generate_report(run_dir: Path | None = None) -> Path:
    """Render the HTML report for a run and write it into the run directory.

    Returns the path to the written ``bioannot_report.html``.
    """
    run_dir = Path(run_dir) if run_dir else resolve_run_dir()
    annotations = _read_annotations(run_dir)
    tsv_order, rows, composite = _read_priority_rows(run_dir)
    genes = _ordered_genes(annotations, tsv_order)
    disease_context = load_disease_context()["context"]
    run_ts = _run_timestamp(run_dir, annotations)

    cyto = _read_cytoscape(run_dir)
    network_payload = _build_network_payload(cyto, annotations, rows, composite)

    html_doc = render_html(
        annotations, genes, rows, composite, network_payload,
        run_dir, disease_context, run_ts,
    )
    out_path = run_dir / "bioannot_report.html"
    out_path.write_text(html_doc, encoding="utf-8")
    return out_path


def main() -> None:
    out = generate_report()
    print(f"Wrote HTML report → {out}")


if __name__ == "__main__":
    main()
