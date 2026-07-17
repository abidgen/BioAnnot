"""NetworkX graph builder and target prioritizer (CLAUDE.md Step 7)."""

from __future__ import annotations

import itertools
import logging
import pickle
from typing import Any

import networkx as nx
import pandas as pd

from src.config import config
from src.utils import load_disease_context

log = logging.getLogger("bio_annot.network")

# The GTEx safety multiplier applied to a flagged target's composite score is
# computed per-tier in src.filters.gtex_safety (config.gtex_tier1_penalty for a
# tier-1 flag, config.gtex_tier2_penalty for tier-2, else 1.00) and carried on
# each record's safety_assessment as ``safety_penalty``; network.py reads that
# value rather than deriving its own.

_NODE_ATTRS = (
    "functions",
    "cellular_states",
    "pathways",
    "disease_associations",
    "druggability_notes",
    "confidence",
    "source_count",
)


def build_target_network(
    final_annotations: dict[str, Any], include_interactor_nodes: bool = True
) -> nx.MultiDiGraph:
    """Build a gene graph from merged annotations.

    A MultiDiGraph so pathway co-membership, direct-interaction, and STRING-PPI
    edges can coexist between the same node pair without overwriting each
    other's attributes. The edge key is the edge type, so each pair carries at
    most one edge of each type per direction.

    Target genes carry ``node_type="target"`` and the full annotation. STRING
    interaction partners (``string_interactors`` on each annotation) add
    ``string_interaction`` edges:
      - between two target genes when STRING links them, and
      - to "satellite" partner nodes (``node_type="interactor"``, no annotation)
        when ``include_interactor_nodes`` is True. Satellites give otherwise
        isolated targets connectivity and let two targets bridge through a
        shared partner; they are excluded from prioritization scoring.
    """
    G = nx.MultiDiGraph()
    targets = set(final_annotations)

    for gene, ann in final_annotations.items():
        G.add_node(
            gene,
            node_type="target",
            functions=ann.get("functions", []),
            cellular_states=ann.get("cellular_states", []),
            pathways=ann.get("pathways", []),
            disease_associations=ann.get("disease_associations", []),
            druggability_notes=ann.get("druggability_notes", ""),
            confidence=ann.get("confidence", 0.0),
            source_count=ann.get("source_count", 0),
            string_interactors=ann.get("string_interactors", []),
            safety_assessment=ann.get("safety_assessment", {}),
            cellxgene_expression=ann.get("cellxgene_expression", {}),
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

    # STRING PPI: weighted by combined score (0–1000). Target↔target edges are
    # deduped to one per unordered pair (STRING partner lists aren't guaranteed
    # symmetric across both genes' top-N, so dedupe by the pair, not by source
    # ordering). Target→satellite edges are added once per (target, partner).
    seen_target_pairs: set[frozenset] = set()
    for gene, ann in final_annotations.items():
        for record in ann.get("string_interactors") or []:
            partner = (record.get("partner") or "").upper()
            score = record.get("combined_score", 0)
            if not partner or partner == gene:
                continue
            if partner in targets:
                pair = frozenset({gene, partner})
                if pair in seen_target_pairs:
                    continue
                seen_target_pairs.add(pair)
                G.add_edge(
                    gene,
                    partner,
                    key="string_interaction",
                    type="string_interaction",
                    weight=score,
                    combined_score=score,
                )
            elif include_interactor_nodes:
                if partner not in G:
                    G.add_node(partner, node_type="interactor")
                if not G.has_edge(gene, partner, key="string_interaction"):
                    G.add_edge(
                        gene,
                        partner,
                        key="string_interaction",
                        type="string_interaction",
                        weight=score,
                        combined_score=score,
                    )

    n_targets = sum(
        1 for _, d in G.nodes(data=True) if d.get("node_type") == "target"
    )
    log.info(
        "Built network: %d nodes (%d targets, %d interactor satellites), %d edges",
        G.number_of_nodes(),
        n_targets,
        G.number_of_nodes() - n_targets,
        G.number_of_edges(),
    )
    return G


def compute_priority_scores(
    G: nx.Graph, disease_filter: str | None = None
) -> list[dict[str, Any]]:
    """Score and rank nodes for target prioritization.

    Composite =
        (w_betweenness·betweenness + w_degree·degree
         + w_disease·min(disease_score, 1.0) + w_druggability·druggability_bonus
         + w_cellxgene·cellxgene_score)
        × confidence × safety_penalty

    The five weights come from src.config (defaults 0.25 / 0.15 / 0.35 / 0.10 /
    0.15, validated to sum to 1.0). ``safety_penalty`` is the per-tier GTEx
    multiplier carried on the record (0.60 for a tier-1 vital-organ flag, 0.80 for
    a tier-2-only flag, else 1.0). ``cellxgene_score`` rewards targets with measured single-cell expression
    (1.0 for ≥3 cell types, 0.5 for ≥1, 0.0 otherwise). disease_score is capped
    at 1.0 in the composite (the raw, uncapped value is still reported).
    """
    # Centrality reflects undirected connectivity: edge direction is nominal
    # (co-membership is symmetric; interaction direction is not meaningful for
    # reach). Project the MultiDiGraph to a simple undirected graph so paths are
    # not constrained by stored edge direction and parallel edges don't distort
    # degree. nx.Graph(G) collapses both direction and multiplicity.
    UG = nx.Graph(G)
    betweenness = nx.betweenness_centrality(UG, normalized=True)
    degree = nx.degree_centrality(UG)

    # Match a disease if any active disease-context scoring term (or the explicit
    # filter) appears in it. Scoring terms come from DISEASE_CONTEXT/DISEASE_TERMS;
    # an optional disease_filter adds one more term (skipped when None/empty so an
    # empty string can't match every disease).
    match_terms = load_disease_context()["scoring_terms"]
    if disease_filter:
        match_terms = match_terms | {disease_filter.lower()}
    strength_weight = {"strong": 1.0, "moderate": 0.5, "weak": 0.2}

    scores: list[dict] = []
    for node, attrs in G.nodes(data=True):
        # Score target genes only; satellite interactor nodes contribute to
        # centrality (above) but are not themselves prioritization candidates.
        if attrs.get("node_type") != "target":
            continue
        disease_assocs = attrs.get("disease_associations", [])
        disease_score = sum(
            strength_weight.get(d.get("evidence_strength", ""), 0.0)
            for d in disease_assocs
            if any(term in d.get("disease", "").lower() for term in match_terms)
        )

        druggability_bonus = 0.2 if attrs.get("druggability_notes") else 0.0
        confidence = attrs.get("confidence", 0.0)

        # CellxGene single-cell expression breadth: reward targets grounded in
        # measured per-cell-type expression. cell_type_count is the number of
        # cell types (after the min-cells filter) the fetcher returned.
        cellxgene = attrs.get("cellxgene_expression") or {}
        cell_type_count = cellxgene.get("cell_type_count", 0)
        cellxgene_score = (
            1.0 if cell_type_count >= 3 else 0.5 if cell_type_count >= 1 else 0.0
        )

        # GTEx safety penalty: deprioritize (don't eliminate) targets with high
        # normal-tissue expression by scaling the composite by the per-tier
        # multiplier already computed on the safety assessment (0.60 tier-1 vital,
        # 0.80 tier-2, else 1.0).
        safety = attrs.get("safety_assessment") or {}
        safety_flag = bool(safety.get("safety_flag", False))
        tier1_flag = bool(safety.get("tier1_flag", False))
        tier2_flag = bool(safety.get("tier2_flag", False))
        safety_penalty = float(safety.get("safety_penalty", 1.0))

        composite = (
            config.weight_betweenness * betweenness[node]
            + config.weight_degree * degree[node]
            + config.weight_disease * min(disease_score, 1.0)
            + config.weight_druggability * druggability_bonus
            + config.weight_cellxgene * cellxgene_score
        ) * confidence * safety_penalty

        scores.append(
            {
                "gene": node,
                "composite": composite,
                "betweenness": betweenness[node],
                "degree": degree[node],
                "disease_score": disease_score,
                "druggability_bonus": druggability_bonus,
                "cellxgene_score": cellxgene_score,
                "confidence": confidence,
                "safety_flag": safety_flag,
                "tier1_flag": tier1_flag,
                "tier2_flag": tier2_flag,
                "safety_penalty": safety_penalty,
                "safety_penalty_applied": safety_penalty < 1.0,
                "tier1_high_tissues": safety.get("tier1_high_tissues", {}),
                "high_expression_tissues": safety.get("high_expression_tissues", []),
                "max_vital_tpm": safety.get("max_vital_tpm", 0.0),
                "max_tpm": safety.get("max_tpm", 0.0),
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


def save_prioritized_tsv(scores: list[dict[str, Any]], path: str) -> None:
    """Write the ranked target table as TSV, flattening list/dict fields."""
    rows = []
    for s in scores:
        rows.append({k: _flatten(v) for k, v in s.items()})
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    log.info("Wrote %d prioritized targets → %s", len(scores), path)
