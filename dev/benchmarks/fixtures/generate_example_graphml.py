"""
Generate synthetic GraphML files for the 2D network tutorial example.

Run from the repo root:
    python examples/generate_example_graphml.py

Writes:
    examples/data/example_specreboot_network.graphml  (30 nodes)
    examples/data/example_phenorewire_network.graphml (22 nodes)
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import networkx as nx

SEED = 42
random.seed(SEED)


def _make_specreboot_network() -> nx.Graph:
    """
    30-node molecular network:
      - 10 nodes with full annotation
      - 10 nodes with partial annotation
      - 10 nodes unannotated
    Edges represent cosine similarity (0.70–1.0 within clusters, random otherwise).
    """
    G = nx.Graph()
    G.name = "specreboot_example"

    annotations = {
        # Fully annotated (metabolite names)
        "F000": "L-Glutamine",
        "F001": "L-Asparagine",
        "F002": "Citric acid",
        "F003": "Isocitric acid",
        "F004": "Succinic acid",
        "F005": "Fumaric acid",
        "F006": "L-Malic acid",
        "F007": "Oxaloacetic acid",
        "F008": "Alpha-Ketoglutaric acid",
        "F009": "Pyruvic acid",
        # Partially annotated (class-level)
        "F010": "amino acid [partial]",
        "F011": "organic acid [partial]",
        "F012": "lipid [partial]",
        "F013": "nucleotide [partial]",
        "F014": "carbohydrate [partial]",
        "F015": None,  # unannotated
        "F016": None,
        "F017": None,
        "F018": None,
        "F019": None,
        # Remaining: unannotated
        **{f"F{i:03d}": None for i in range(20, 30)},
    }

    for nid, ann in annotations.items():
        attrs = {
            "feature_id": nid,
            "label": ann or nid,
        }
        if ann:
            attrs["consensus_annotation"] = ann
        G.add_node(nid, **attrs)

    rng = random.Random(SEED)

    # Cluster 1: F000-F004 (amino acids / TCA cycle entry)
    cluster1 = ["F000", "F001", "F002", "F003", "F004"]
    for i, u in enumerate(cluster1):
        for v in cluster1[i + 1:]:
            G.add_edge(u, v, cosine=round(rng.uniform(0.78, 0.98), 3), weight=round(rng.uniform(0.78, 0.98), 3))

    # Cluster 2: F005-F009 (TCA cycle intermediates)
    cluster2 = ["F005", "F006", "F007", "F008", "F009"]
    for i, u in enumerate(cluster2):
        for v in cluster2[i + 1:]:
            G.add_edge(u, v, cosine=round(rng.uniform(0.72, 0.95), 3), weight=round(rng.uniform(0.72, 0.95), 3))

    # Bridge between clusters
    G.add_edge("F004", "F005", cosine=0.71, weight=0.71)

    # Partial annotation nodes connected to clusters
    G.add_edge("F010", "F000", cosine=round(rng.uniform(0.70, 0.82), 3), weight=0.75)
    G.add_edge("F011", "F005", cosine=round(rng.uniform(0.70, 0.82), 3), weight=0.73)
    G.add_edge("F012", "F002", cosine=round(rng.uniform(0.70, 0.78), 3), weight=0.72)
    G.add_edge("F013", "F008", cosine=round(rng.uniform(0.70, 0.80), 3), weight=0.74)
    G.add_edge("F014", "F009", cosine=round(rng.uniform(0.70, 0.79), 3), weight=0.71)

    # Unannotated nodes sparsely connected
    for i in range(15, 20):
        nbr = f"F{rng.randint(0, 14):03d}"
        G.add_edge(f"F{i:03d}", nbr, cosine=round(rng.uniform(0.70, 0.80), 3), weight=0.72)

    for i in range(20, 30):
        nbr = f"F{rng.randint(0, 19):03d}"
        G.add_edge(f"F{i:03d}", nbr, cosine=round(rng.uniform(0.70, 0.75), 3), weight=0.71)

    return G


def _make_phenorewire_network() -> nx.Graph:
    """
    22-node PhenoRewire rewiring network (same feature IDs as example dataset).
    Nodes include rewiring_score and triage_score attributes.
    """
    G = nx.Graph()
    G.name = "phenorewire_example"

    rng = random.Random(SEED + 1)

    feature_ids = [f"F{i:03d}" for i in range(22)]
    for fid in feature_ids:
        G.add_node(
            fid,
            feature_id=fid,
            rewiring_score=round(rng.uniform(0.0, 1.0), 3),
            triage_score=round(rng.uniform(0.3, 0.95), 3),
            state="rewired" if rng.random() > 0.4 else "shared",
        )

    # Rewired edges (REF-only and CASE-only)
    ref_only = [("F000", "F001"), ("F002", "F003"), ("F004", "F005"), ("F006", "F007")]
    case_only = [("F001", "F002"), ("F003", "F004"), ("F005", "F006"), ("F008", "F009")]
    shared = [("F000", "F009"), ("F010", "F011"), ("F012", "F013"), ("F014", "F015")]
    sign_switches = [("F016", "F017"), ("F018", "F019")]

    for u, v in ref_only:
        G.add_edge(u, v, state="ref_only", r=round(rng.uniform(0.55, 0.90), 3))
    for u, v in case_only:
        G.add_edge(u, v, state="case_only", r=round(rng.uniform(0.55, 0.90), 3))
    for u, v in shared:
        G.add_edge(u, v, state="shared", r=round(rng.uniform(0.60, 0.95), 3))
    for u, v in sign_switches:
        G.add_edge(u, v, state="sign_switch", r=round(rng.uniform(0.55, 0.85), 3))

    return G


def main() -> None:
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    specreboot = _make_specreboot_network()
    pr_path = out_dir / "example_specreboot_network.graphml"
    nx.write_graphml(specreboot, str(pr_path))
    print(f"Written: {pr_path}  ({specreboot.number_of_nodes()} nodes, {specreboot.number_of_edges()} edges)")

    phenorewire = _make_phenorewire_network()
    mn_path = out_dir / "example_phenorewire_network.graphml"
    nx.write_graphml(phenorewire, str(mn_path))
    print(f"Written: {mn_path}  ({phenorewire.number_of_nodes()} nodes, {phenorewire.number_of_edges()} edges)")


if __name__ == "__main__":
    main()
