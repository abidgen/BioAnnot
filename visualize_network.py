"""Visualize the target network and prioritization scores.

Standalone — reads existing pipeline outputs only, no pipeline rerun:
  - outputs/target_network.gpickle   (NetworkX MultiDiGraph)
  - outputs/prioritized_targets.tsv  (composite scores + components)

Generates four plots into outputs/plots/:
  a) target_network.png        — graph, nodes by composite score, edges by type
  b) score_breakdown.png       — per-gene confidence / disease / network / composite
  c) pathway_heatmap.png       — genes × canonical pathways presence matrix
  d) cellxgene_expression.png  — top 5 CellxGene cell types per gene (if present)

Run with: python visualize_network.py
"""

from __future__ import annotations

import pickle
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display needed

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

NETWORK_PATH = "outputs/target_network.gpickle"
SCORES_PATH = "outputs/prioritized_targets.tsv"
PLOTS_DIR = Path("outputs/plots")

NON_CANONICAL_PREFIX = "NON-CANONICAL: "

EDGE_COLORS = {
    "pathway_comembership": "green",
    "direct_interaction": "orange",
    "string_interaction": "steelblue",
}

# Per-edge-type curvature: parallel edges between the same node pair (e.g. a
# direct_interaction and a string_interaction between TP53 and BRCA1) are drawn
# with different curvature so they don't overlap and both stay visible.
EDGE_CURVATURE = {
    "pathway_comembership": "arc3,rad=0.0",
    "direct_interaction": "arc3,rad=-0.2",
    "string_interaction": "arc3,rad=0.2",
}

# Node-size range (points^2) mapped across the composite-score range.
NODE_SIZE_MIN = 600
NODE_SIZE_MAX = 3000

# Deterministic spring-layout seed (shown in the plot subtitle) and the minimum
# spacing enforced between nodes after layout to prevent overlap.
LAYOUT_SEED = 42
MIN_NODE_DIST = 0.3


def load_network(path: str) -> nx.MultiDiGraph:
    """Load the pickled NetworkX graph."""
    with open(path, "rb") as f:
        return pickle.load(f)


def load_scores(path: str) -> pd.DataFrame:
    """Load the prioritized-targets TSV indexed by gene symbol."""
    df = pd.read_csv(path, sep="\t")
    return df.set_index("gene")


def _composite_map(scores: pd.DataFrame) -> dict[str, float]:
    """Map gene -> composite score (0.0 for any node missing from the table)."""
    return scores["composite"].to_dict()


def target_subgraph(G: nx.MultiDiGraph) -> nx.MultiDiGraph:
    """Subgraph of target nodes and the edges among them.

    STRING enrichment can add many satellite interactor nodes (node_type
    "interactor"); plotting them produces an unreadable hairball, so the
    visualizations focus on the target genes and their mutual relationships.
    Falls back to the whole graph if no node_type is set (older pickles).
    """
    targets = [n for n, d in G.nodes(data=True) if d.get("node_type") == "target"]
    if not targets:
        return G
    return G.subgraph(targets).copy()


def _enforce_min_distance(
    pos: dict, min_dist: float = MIN_NODE_DIST, iterations: int = 100
) -> dict:
    """Push apart any node pair closer than ``min_dist`` (deterministic).

    Spring layout can leave nodes overlapping on small graphs; this is a simple
    relaxation pass that nudges too-close pairs apart along their connecting axis
    and repeats until everything is separated (or ``iterations`` is hit). No
    randomness, so the result stays reproducible for a given input layout.
    """
    keys = list(pos)
    points = {k: np.asarray(v, dtype=float) for k, v in pos.items()}
    for _ in range(iterations):
        moved = False
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                a, b = keys[i], keys[j]
                delta = points[a] - points[b]
                dist = float(np.linalg.norm(delta))
                if dist >= min_dist:
                    continue
                if dist < 1e-9:
                    # Coincident nodes: separate along a fixed axis deterministically.
                    delta = np.array([1.0, 0.0])
                    dist = 1.0
                shift = (min_dist - dist) / 2.0
                direction = delta / dist
                points[a] += direction * shift
                points[b] -= direction * shift
                moved = True
        if not moved:
            break
    return {k: tuple(v) for k, v in points.items()}


def plot_target_network(
    G: nx.MultiDiGraph, scores: pd.DataFrame, out_path: Path
) -> None:
    """Draw the gene graph: node color/size by composite, edges colored by type."""
    G = target_subgraph(G)  # exclude STRING satellite nodes from the plot
    composite = _composite_map(scores)
    nodes = list(G.nodes())
    values = np.array([composite.get(n, 0.0) for n in nodes])

    vmin, vmax = float(values.min()), float(values.max())
    norm = Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1e-9)
    cmap = plt.get_cmap("coolwarm")  # blue (low) -> red (high)
    node_colors = [cmap(norm(v)) for v in values]

    # Size proportional to composite score (linear across the observed range).
    if vmax > vmin:
        sizes = [
            NODE_SIZE_MIN + (NODE_SIZE_MAX - NODE_SIZE_MIN) * (v - vmin) / (vmax - vmin)
            for v in values
        ]
    else:
        sizes = [(NODE_SIZE_MIN + NODE_SIZE_MAX) / 2] * len(nodes)

    pos = nx.spring_layout(G, seed=LAYOUT_SEED, k=2.5)
    pos = _enforce_min_distance(pos, min_dist=MIN_NODE_DIST)

    fig, ax = plt.subplots(figsize=(14, 10))

    nx.draw_networkx_nodes(
        G, pos, nodelist=nodes, node_color=node_colors, node_size=sizes, ax=ax
    )
    nx.draw_networkx_labels(G, pos, font_size=11, font_weight="bold", ax=ax)

    # Draw edges grouped by type so each gets its own color. Edge type is conveyed
    # by color via the legend below; per-edge text labels are omitted as redundant.
    for etype, color in EDGE_COLORS.items():
        edges = [
            (u, v) for u, v, d in G.edges(data=True) if d.get("type") == etype
        ]
        if not edges:
            continue
        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=edges,
            edge_color=color,
            width=2.0,
            alpha=0.8,
            connectionstyle=EDGE_CURVATURE.get(etype, "arc3,rad=0.1"),
            ax=ax,
        )

    # Legend for edge types.
    legend_handles = [
        plt.Line2D([0], [0], color=color, lw=2, label=etype)
        for etype, color in EDGE_COLORS.items()
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        title="Edge type",
    )

    # Colorbar for composite score.
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label("composite score")

    ax.set_title(
        "Target Network — node color/size by composite score\n"
        f"(layout seed={LAYOUT_SEED})"
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_score_breakdown(scores: pd.DataFrame, out_path: Path) -> None:
    """Grouped horizontal bars per gene: confidence, disease, network, composite.

    network_score is the sum of the centrality components (betweenness + degree),
    since the TSV has no single network column.
    """
    df = scores.copy()
    df["network_score"] = df["betweenness"] + df["degree"]
    # disease_score is raw and unbounded (can exceed 6), but the composite uses
    # min(disease_score, 1.0). Plot the capped contribution so all four metrics
    # share a 0-1 scale and the bars are comparable.
    df["disease_capped"] = df["disease_score"].clip(upper=1.0)

    metrics = ["confidence", "disease_capped", "network_score", "composite"]
    labels = [
        "confidence",
        "disease_score (capped@1.0)",
        "network_score",
        "composite_score",
    ]

    # Order genes by composite ascending so the top target sits at the top.
    df = df.sort_values("composite", ascending=True)
    genes = list(df.index)

    y = np.arange(len(genes))
    bar_h = 0.8 / len(metrics)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.9 * len(genes) + 2)))
    cmap = plt.get_cmap("tab10")
    for i, (metric, label) in enumerate(zip(metrics, labels)):
        offset = (i - (len(metrics) - 1) / 2) * bar_h
        ax.barh(
            y + offset,
            df[metric].values,
            height=bar_h,
            label=label,
            color=cmap(i),
        )

    ax.set_yticks(y)
    ax.set_yticklabels(genes)
    ax.set_xlabel("score")
    ax.set_title("Per-gene score breakdown")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_pathway_heatmap(G: nx.MultiDiGraph, out_path: Path) -> None:
    """Heatmap of genes × canonical pathways (1 = gene has pathway)."""
    G = target_subgraph(G)  # satellite interactor nodes carry no pathways
    gene_pathways: dict[str, set[str]] = {}
    for gene, attrs in G.nodes(data=True):
        canonical = {
            p
            for p in attrs.get("pathways", [])
            if not p.startswith(NON_CANONICAL_PREFIX)
        }
        gene_pathways[gene] = canonical

    all_pathways = sorted(set().union(*gene_pathways.values())) if gene_pathways else []
    genes = sorted(gene_pathways)

    if not all_pathways:
        # Nothing canonical to show — emit a placeholder so the file still exists.
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No canonical pathways found", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    matrix = np.array(
        [[1 if p in gene_pathways[g] else 0 for p in all_pathways] for g in genes]
    )

    # Truncate long Reactome names so rotated tick labels stay legible.
    def _short(name: str, limit: int = 45) -> str:
        return name if len(name) <= limit else name[: limit - 1] + "…"

    # Width scales with pathway count; height with gene count.
    fig_w = max(8, 0.30 * len(all_pathways) + 4)
    fig_h = max(4, 0.5 * len(genes) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(all_pathways)))
    ax.set_xticklabels([_short(p) for p in all_pathways], rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes)
    ax.set_title("Canonical pathway membership (genes × pathways)")

    # Thin grid between cells for readability.
    ax.set_xticks(np.arange(-0.5, len(all_pathways), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="lightgray", linewidth=0.5)
    ax.tick_params(which="minor", length=0)

    # Reserve generous bottom space for the long rotated labels (tight_layout
    # can't always fit them), and avoid the UserWarning it raises when it can't.
    fig.subplots_adjust(bottom=0.55, left=0.08, right=0.98, top=0.93)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_cellxgene_expression(G: nx.MultiDiGraph, out_path: Path) -> None:
    """Top-5 CellxGene cell types per gene, one horizontal-bar subplot per gene.

    Reads the ``cellxgene_expression`` node attribute (written by the pipeline);
    only genes that carry it with a non-empty ``top_cell_types`` are plotted.
    """
    genes = [
        (n, d["cellxgene_expression"])
        for n, d in G.nodes(data=True)
        if d.get("node_type") == "target"
        and (d.get("cellxgene_expression") or {}).get("top_cell_types")
    ]

    if not genes:
        # Nothing to show (CellxGene disabled or no data) — emit a placeholder.
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No CellxGene expression data", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    def _short(name: str, limit: int = 40) -> str:
        return name if len(name) <= limit else name[: limit - 1] + "…"

    n = len(genes)
    fig, axes = plt.subplots(n, 1, figsize=(9, max(3, 2.6 * n)), squeeze=False)
    cmap = plt.get_cmap("viridis")
    for ax, (gene, cx) in zip(axes[:, 0], genes):
        top = cx["top_cell_types"][:5]
        cell_types = [t["cell_type"] for t in top]
        means = [t["mean_expr"] for t in top]
        # Highest-expression cell type at the top of the subplot.
        ypos = np.arange(len(cell_types))[::-1]
        colors = cmap(np.linspace(0.25, 0.85, len(cell_types)))
        ax.barh(ypos, means, color=colors)
        ax.set_yticks(ypos)
        ax.set_yticklabels([_short(c) for c in cell_types], fontsize=8)
        ax.set_xlabel("mean expression")
        ax.set_title(
            f"{gene} — top cell types ({cx.get('tissue', '')})", fontsize=10
        )

    fig.suptitle("CellxGene Census — top 5 cell types per gene")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def print_summary(G: nx.MultiDiGraph, scores: pd.DataFrame) -> None:
    """Print node/edge counts and the top gene by composite score."""
    edge_types = Counter(d.get("type", "unknown") for _, _, d in G.edges(data=True))

    print("=== Network summary ===")
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")
    for etype, count in sorted(edge_types.items()):
        print(f"  {etype}: {count}")

    top_gene = scores["composite"].idxmax()
    top_score = scores["composite"].max()
    print(f"Top gene by composite score: {top_gene} ({top_score:.4f})")


def main() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    G = load_network(NETWORK_PATH)
    scores = load_scores(SCORES_PATH)

    plot_target_network(G, scores, PLOTS_DIR / "target_network.png")
    plot_score_breakdown(scores, PLOTS_DIR / "score_breakdown.png")
    plot_pathway_heatmap(G, PLOTS_DIR / "pathway_heatmap.png")
    plot_cellxgene_expression(G, PLOTS_DIR / "cellxgene_expression.png")

    print(f"Saved plots to {PLOTS_DIR}/")
    print("  - target_network.png")
    print("  - score_breakdown.png")
    print("  - pathway_heatmap.png")
    print("  - cellxgene_expression.png")
    print()
    print_summary(G, scores)

    # Also export the network to Cytoscape.js JSON (full + targets-only)
    # whenever plots are generated.
    from scripts.export_cytoscape import export_cytoscape

    export_cytoscape(G)


if __name__ == "__main__":
    main()
