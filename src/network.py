"""NetworkX graph builder and target prioritizer (CLAUDE.md Step 7)."""

from __future__ import annotations

import itertools
import logging
import pickle

import networkx as nx
import pandas as pd

log = logging.getLogger("bio_annot.network")

# Oncology term set used to match a disease against the disease_filter.
# A disease counts toward disease_score if any of these terms appears in it.
ONCOLOGY_TERMS = {
    "cancer",
    "carcinoma",
    "sarcoma",
    "lymphoma",
    "leukemia",
    "melanoma",
    "adenocarcinoma",
    "glioma",
    "myeloma",
    "blastoma",
    "tumor",
    "tumour",
    "neoplasm",
    "malignancy",
    "oncology",
}

_NODE_ATTRS = (
    "functions",
    "cellular_states",
    "pathways",
    "disease_associations",
    "druggability_notes",
    "confidence",
    "source_count",
)


def build_target_network(final_annotations: dict) -> nx.MultiDiGraph:
    """Build a gene graph from merged annotations.

    A MultiDiGraph so pathway co-membership and direct-interaction edges can
    coexist between the same node pair without overwriting each other's
    attributes. The edge key is the edge type, so each pair carries at most one
    edge of each type per direction.
    """
    G = nx.MultiDiGraph()

    for gene, ann in final_annotations.items():
        G.add_node(
            gene,
            functions=ann.get("functions", []),
            cellular_states=ann.get("cellular_states", []),
            pathways=ann.get("pathways", []),
            disease_associations=ann.get("disease_associations", []),
            druggability_notes=ann.get("druggability_notes", ""),
            confidence=ann.get("confidence", 0.0),
            source_count=ann.get("source_count", 0),
        )

    # Pathway co-membership: any two genes sharing ≥1 pathway. Symmetric, so
    # stored as one directed edge per unordered pair (direction is nominal).
    pathway_sets = {
        gene: set(ann.get("pathways", []))
        for gene, ann in final_annotations.items()
    }
    for gene_a, gene_b in itertools.combinations(final_annotations, 2):
        shared = pathway_sets[gene_a] & pathway_sets[gene_b]
        if shared:
            shared_list = sorted(shared)
            G.add_edge(
                gene_a,
                gene_b,
                key="pathway_comembership",
                type="pathway_comembership",
                shared_pathways=shared_list,
                weight=len(shared_list),
            )

    # Direct interaction: directed edge gene → any interactor that is a node.
    for gene, ann in final_annotations.items():
        for interactor in ann.get("interactors") or []:
            if interactor == gene or interactor not in G:
                continue
            G.add_edge(
                gene,
                interactor,
                key="direct_interaction",
                type="direct_interaction",
                weight=2,
                source=gene,
            )

    log.info(
        "Built network: %d nodes, %d edges", G.number_of_nodes(), G.number_of_edges()
    )
    return G


def compute_priority_scores(G: nx.Graph, disease_filter: str) -> list[dict]:
    """Score and rank nodes for target prioritization.

    Composite =
        (0.30·betweenness + 0.20·degree + 0.40·min(disease_score, 1.0)
         + 0.10·druggability_bonus) × confidence
    """
    betweenness = nx.betweenness_centrality(G, normalized=True)
    degree = nx.degree_centrality(G)

    # Match a disease if any oncology term (or the explicit filter) appears in it.
    match_terms = ONCOLOGY_TERMS | {disease_filter.lower()}
    strength_weight = {"strong": 1.0, "moderate": 0.5, "weak": 0.2}

    scores: list[dict] = []
    for node, attrs in G.nodes(data=True):
        disease_assocs = attrs.get("disease_associations", [])
        disease_score = sum(
            strength_weight.get(d.get("evidence_strength", ""), 0.0)
            for d in disease_assocs
            if any(term in d.get("disease", "").lower() for term in match_terms)
        )

        druggability_bonus = 0.2 if attrs.get("druggability_notes") else 0.0
        confidence = attrs.get("confidence", 0.0)

        composite = (
            0.30 * betweenness[node]
            + 0.20 * degree[node]
            + 0.40 * min(disease_score, 1.0)
            + 0.10 * druggability_bonus
        ) * confidence

        scores.append(
            {
                "gene": node,
                "composite": composite,
                "betweenness": betweenness[node],
                "degree": degree[node],
                "disease_score": disease_score,
                "druggability_bonus": druggability_bonus,
                "confidence": confidence,
                "pathways": attrs.get("pathways", []),
                "disease_associations": disease_assocs,
            }
        )

    scores.sort(key=lambda s: s["composite"], reverse=True)
    return scores


def save_network(G: nx.Graph, path: str) -> None:
    """Persist the graph to a gpickle file.

    NetworkX removed nx.write_gpickle in 3.0 (this project pins networkx>=3.3),
    so fall back to the stdlib pickle, which is what gpickle always used.
    """
    if hasattr(nx, "write_gpickle"):
        nx.write_gpickle(G, path)
    else:
        with open(path, "wb") as f:
            pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    log.info("Wrote network → %s", path)


def _flatten(value) -> str:
    """Flatten a list/dict field into a pipe-separated string for TSV output."""
    if isinstance(value, list):
        parts = []
        for item in value:
            if isinstance(item, dict):
                parts.append(
                    ":".join(str(v) for v in item.values())
                )
            else:
                parts.append(str(item))
        return "|".join(parts)
    return str(value)


def save_prioritized_tsv(scores: list[dict], path: str) -> None:
    """Write the ranked target table as TSV, flattening list/dict fields."""
    rows = []
    for s in scores:
        rows.append({k: _flatten(v) for k, v in s.items()})
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    log.info("Wrote %d prioritized targets → %s", len(scores), path)
