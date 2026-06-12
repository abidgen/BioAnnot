"""Visualize the target network and prioritization scores.

Standalone — reads existing pipeline outputs only, no pipeline rerun:
  - outputs/target_network.gpickle   (NetworkX MultiDiGraph)
  - outputs/prioritized_targets.tsv  (composite scores + components)

Generates three plots into outputs/plots/:
  a) target_network.png   — graph, nodes by composite score, edges by type
  b) score_breakdown.png  — per-gene confidence / disease / network / composite
  c) pathway_heatmap.png  — genes × canonical pathways presence matrix

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
}

# Node-size range (points^2) mapped across the composite-score range.
NODE_SIZE_MIN = 600
NODE_SIZE_MAX = 3000


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


def plot_target_network(
    G: nx.MultiDiGraph, scores: pd.DataFrame, out_path: Path
) -> None:
    """Draw the gene graph: node color/size by composite, edges colored by type."""
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

    pos = nx.spring_layout(G, seed=42, k=1.5)

    fig, ax = plt.subplots(figsize=(11, 9))

    nx.draw_networkx_nodes(
        G, pos, nodelist=nodes, node_color=node_colors, node_size=sizes, ax=ax
    )
    nx.draw_networkx_labels(G, pos, font_size=11, font_weight="bold", ax=ax)

    # Draw edges grouped by type so each gets its own color, and collect labels.
    edge_labels: dict[tuple[str, str], str] = {}
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
            connectionstyle="arc3,rad=0.08",
            ax=ax,
        )
        for u, v in edges:
            edge_labels[(u, v)] = etype

    nx.draw_networkx_edge_labels(
        G, pos, edge_labels=edge_labels, font_size=7, font_color="dimgray", ax=ax
    )

    # Legend for edge types.
    legend_handles = [
        plt.Line2D([0], [0], color=color, lw=2, label=etype)
        for etype, color in EDGE_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="lower left", title="Edge type")

    # Colorbar for composite score.
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label("composite score")

    ax.set_title("Target Network — node color/size by composite score")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
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

    print(f"Saved plots to {PLOTS_DIR}/")
    print("  - target_network.png")
    print("  - score_breakdown.png")
    print("  - pathway_heatmap.png")
    print()
    print_summary(G, scores)


if __name__ == "__main__":
    main()
