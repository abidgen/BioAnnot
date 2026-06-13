"""Export the target network to Cytoscape.js JSON (standalone, post-processing).

Produces two files from a network graph via nx.cytoscape_data(G):
  - outputs/target_network_cytoscape.json          full graph (all nodes +
    STRING satellite interactors) — for deep PPI exploration
  - outputs/target_network_cytoscape_targets.json  target nodes + their mutual
    edges only — for clean presentation

Both Cytoscape and Cytoscape.js can import these directly; node/edge attributes
are carried through under each element's "data". Run standalone (loads
outputs/target_network.gpickle), or call export_cytoscape(G) from
visualize_network.py:

    python scripts/export_cytoscape.py
"""

from __future__ import annotations

import json
import pickle
from collections import Counter
from pathlib import Path

import networkx as nx

NETWORK_PATH = Path("outputs/target_network.gpickle")
OUTPUT_PATH = Path("outputs/target_network_cytoscape.json")
TARGETS_OUTPUT_PATH = Path("outputs/target_network_cytoscape_targets.json")
FULL_CX2_PATH = Path("outputs/target_network_cytoscape.cx2")
TARGETS_CX2_PATH = Path("outputs/target_network_cytoscape_targets.cx2")

# Cytoscape's own per-element keys (added by nx.cytoscape_data), excluded from
# the "domain attributes" line of the summary.
_CYTOSCAPE_KEYS = {"id", "value", "name", "source", "target", "key"}


def _write(data: dict, path: Path) -> None:
    """Write a cytoscape_data dict to JSON (default=str guards odd values)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)


def _summarize(label: str, data: dict, path: Path) -> None:
    """Print node/edge counts (by type) for one exported file."""
    nodes = data["elements"].get("nodes", [])
    edges = data["elements"].get("edges", [])
    node_types = Counter(n.get("data", {}).get("node_type", "?") for n in nodes)
    edge_types = Counter(e.get("data", {}).get("type", "?") for e in edges)
    print(f"  {label} → {path}")
    print(f"      nodes: {len(nodes)}  by node_type: {dict(node_types)}")
    print(f"      edges: {len(edges)}  by type: {dict(edge_types)}")


def _cx2_scalar_ok(value) -> bool:
    """True if a value is a CX2-native scalar or a list of scalars.

    CX2 attribute values must be scalars or lists of scalars; nested dicts and
    lists-of-dicts (e.g. cellxgene_expression, disease_associations) are not
    allowed and must be serialized to a JSON string first.
    """
    if value is None or isinstance(value, (str, int, float, bool)):
        return True
    if isinstance(value, list):
        return all(isinstance(x, (str, int, float, bool)) for x in value)
    return False


def _sanitize_for_cx2(G: nx.Graph) -> nx.Graph:
    """Return a copy of G with non-CX2-native node/edge attrs JSON-stringified."""
    H = G.copy()
    for _node, attrs in H.nodes(data=True):
        for key, value in list(attrs.items()):
            if not _cx2_scalar_ok(value):
                attrs[key] = json.dumps(value, default=str)
    for _u, _v, attrs in H.edges(data=True):
        for key, value in list(attrs.items()):
            if not _cx2_scalar_ok(value):
                attrs[key] = json.dumps(value, default=str)
    return H


def _export_cx2(graph: nx.Graph, path: Path, label: str) -> bool:
    """Export a graph to a CX2 file via ndex2; skip if ndex2 is unavailable.

    ndex2 3.x has no NiceCXNetwork.export_to_cx2; the supported path is the CX2
    factory (NetworkXToCX2NetworkFactory → CX2Network.write_as_raw_cx2).
    """
    try:
        from ndex2.cx2 import NetworkXToCX2NetworkFactory
    except ImportError:
        print(f"  cx2: ndex2 not installed; skipping {path} (pip install ndex2)")
        return False

    cx2net = NetworkXToCX2NetworkFactory().get_cx2network(_sanitize_for_cx2(graph))
    path.parent.mkdir(parents=True, exist_ok=True)
    cx2net.write_as_raw_cx2(str(path))
    print(
        f"  cx2 ({label}) → {path}  "
        f"({len(cx2net.get_nodes())} nodes, {len(cx2net.get_edges())} edges)"
    )
    return True


def export_cytoscape(
    G: nx.Graph,
    output_path: Path = OUTPUT_PATH,
    targets_only_path: Path = TARGETS_OUTPUT_PATH,
    full_cx2_path: Path = FULL_CX2_PATH,
    targets_cx2_path: Path = TARGETS_CX2_PATH,
) -> tuple[dict, dict]:
    """Export ``G`` to Cytoscape.js JSON (full + targets-only) and matching CX2.

    Returns ``(full_data, targets_data)``.
    """
    output_path = Path(output_path)
    targets_only_path = Path(targets_only_path)
    full_cx2_path = Path(full_cx2_path)
    targets_cx2_path = Path(targets_cx2_path)

    # Full graph export (all nodes, including STRING satellite interactors).
    data = nx.cytoscape_data(G)
    _write(data, output_path)

    # Target-only export: target nodes and the edges among them.
    target_nodes = [
        n for n, d in G.nodes(data=True) if d.get("node_type") == "target"
    ]
    G_targets = G.subgraph(target_nodes).copy()
    data_targets = nx.cytoscape_data(G_targets)
    _write(data_targets, targets_only_path)

    domain_attrs = set()
    for node in data["elements"].get("nodes", []):
        domain_attrs.update(node.get("data", {}))
    domain_attrs -= _CYTOSCAPE_KEYS

    print("Exported Cytoscape JSON (2 files):")
    _summarize("full   ", data, output_path)
    _summarize("targets", data_targets, targets_only_path)
    print(f"  node attributes included: {sorted(domain_attrs)}")

    # CX2 exports (for NDEx / Cytoscape CX2 import): full graph and targets only.
    _export_cx2(G, full_cx2_path, "full")
    _export_cx2(G_targets, targets_cx2_path, "targets")
    return data, data_targets


def _load_network(path: Path = NETWORK_PATH) -> nx.Graph | None:
    """Load the pickled network graph, or None if the file is missing."""
    path = Path(path)
    if not path.exists():
        print(f"No network at {path}; run the pipeline first.")
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


if __name__ == "__main__":
    graph = _load_network()
    if graph is not None:
        export_cytoscape(graph)
