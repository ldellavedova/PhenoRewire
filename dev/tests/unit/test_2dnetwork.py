"""
Unit tests for phenorewire.twodn.chemical_integration.

Tests:
  - 1-hop and 2-hop neighborhood expansion
  - min_cosine filter correctly excludes low-similarity edges
  - singleton handling (feature not in MN → zero neighbors)
  - feature not in MN at all → handled gracefully
  - annotation_confidence computation
  - all required output columns present in both tables
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

from phenorewire.twodn.chemical_integration import build_annotation_expansion


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_graphml(G: nx.Graph, path: Path) -> Path:
    nx.write_graphml(G, str(path))
    return path


def _make_mn_graph() -> nx.Graph:
    """
    Chemical similarity network:

        A --0.90-- B --0.80-- C --0.75-- D
                   |
                  0.65
                   |
                   E --0.72-- F(annotated)

    Annotations: B="CompoundB", F="CompoundF", D="CompoundD"
    """
    G = nx.Graph()
    for n, ann in [("A", None), ("B", "CompoundB"), ("C", None), ("D", "CompoundD"),
                   ("E", None), ("F", "CompoundF")]:
        attrs = {"feature_id": n}
        if ann:
            attrs["consensus_annotation"] = ann
        G.add_node(n, **attrs)
    G.add_edge("A", "B", cosine=0.90, weight=0.90)
    G.add_edge("B", "C", cosine=0.80, weight=0.80)
    G.add_edge("C", "D", cosine=0.75, weight=0.75)
    G.add_edge("B", "E", cosine=0.65, weight=0.65)
    G.add_edge("E", "F", cosine=0.72, weight=0.72)
    return G


def _make_pr_graph(node_ids: list[str]) -> nx.Graph:
    """Minimal PhenoRewire network for the same feature IDs."""
    G = nx.Graph()
    for n in node_ids:
        G.add_node(n, feature_id=n, rewiring_score=0.5, triage_score=0.6)
    for u, v in zip(node_ids, node_ids[1:]):
        G.add_edge(u, v, state="shared")
    return G


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def graphml_files(tmp_path_factory):
    base = tmp_path_factory.mktemp("2dnetwork_test")
    mn = _make_mn_graph()
    pr = _make_pr_graph(["A", "B", "C"])
    mn_path = base / "mn.graphml"
    pr_path = base / "pr.graphml"
    _write_graphml(mn, mn_path)
    _write_graphml(pr, pr_path)
    return mn_path, pr_path, base


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_1hop_neighbors_correct(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["A"]})
    out = base / "out_1hop"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=1,
        min_cosine=0.70,
        output_dir=out,
    )
    assert not summary.empty
    row = summary[summary["feature_id"] == "A"].iloc[0]
    # A has one 1-hop neighbor: B (cosine=0.90 >= 0.70)
    assert row["n_mn_neighbors"] == 1
    assert row["n_annotated_neighbors"] == 1  # B is annotated
    assert float(row["annotation_confidence"]) == pytest.approx(1.0)


def test_2hop_neighbors_correct(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["A"]})
    out = base / "out_2hop"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=2,
        min_cosine=0.70,
        output_dir=out,
    )
    row = summary[summary["feature_id"] == "A"].iloc[0]
    # 2-hop from A: B (hop1, cos=0.90) and C (hop2, cos=0.80)
    # E is hop2 but via B→E(0.65) which is below min_cosine=0.70 → excluded
    assert row["n_mn_neighbors"] == 2
    neighbor_ids = expansion[expansion["feature_id"] == "A"]["neighbor_id"].tolist()
    assert "B" in neighbor_ids
    assert "C" in neighbor_ids
    assert "E" not in neighbor_ids


def test_min_cosine_filter(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["B"]})
    out = base / "out_mincosine"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=2,
        min_cosine=0.70,  # B->E edge (0.65) should be excluded
        output_dir=out,
    )
    row = summary[summary["feature_id"] == "B"].iloc[0]
    # From B with min_cosine=0.70:
    #   hop1: A(0.90 ok), C(0.80 ok), E(0.65 excluded)
    #   hop2 from A: none (A has only B as neighbor)
    #   hop2 from C: D(0.75 ok)
    # Total: A, C, D = 3 neighbors
    assert row["n_mn_neighbors"] == 3
    neighbor_ids = expansion[expansion["feature_id"] == "B"]["neighbor_id"].tolist()
    assert "E" not in neighbor_ids
    assert "A" in neighbor_ids
    assert "D" in neighbor_ids


def test_singleton_not_in_mn(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["NOTINMN"]})
    out = base / "out_singleton"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=2,
        min_cosine=0.70,
        output_dir=out,
    )
    assert not summary.empty
    row = summary[summary["feature_id"] == "NOTINMN"].iloc[0]
    assert row["n_mn_neighbors"] == 0
    assert row["in_specreboot"] is False or row["in_specreboot"] == False
    assert row["annotation_confidence"] == 0.0


def test_annotation_confidence_mixed(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["C"]})
    out = base / "out_conf"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=1,
        min_cosine=0.70,
        output_dir=out,
    )
    row = summary[summary["feature_id"] == "C"].iloc[0]
    # From C, hop1: B (annotated, cos=0.80), D (annotated, cos=0.75)
    # annotation_confidence = 2/2 = 1.0
    assert float(row["annotation_confidence"]) == pytest.approx(1.0)
    assert row["n_annotated_neighbors"] == 2


def test_output_files_exist(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["A", "B"]})
    out = base / "out_files"
    build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=2,
        min_cosine=0.70,
        output_dir=out,
    )
    assert (out / "annotation_summary.csv").exists()
    assert (out / "annotation_expansion.csv").exists()
    assert (out / "annotation_first_neighbours.csv").exists()
    assert (out / "subnetworks").is_dir()


def test_summary_required_columns(graphml_files):
    from phenorewire.twodn.chemical_integration import _SUMMARY_COLS
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["A"]})
    out = base / "out_cols"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=1,
        min_cosine=0.70,
        output_dir=out,
    )
    for col in _SUMMARY_COLS:
        assert col in summary.columns, f"Missing summary column: {col}"


def test_expansion_required_columns(graphml_files):
    from phenorewire.twodn.chemical_integration import _EXPANSION_COLS
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["A"]})
    out = base / "out_exp_cols"
    summary, expansion, _ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=pd.DataFrame(),
        hop_depth=2,
        min_cosine=0.70,
        output_dir=out,
    )
    if not expansion.empty:
        for col in _EXPANSION_COLS:
            assert col in expansion.columns, f"Missing expansion column: {col}"


def test_rewiring_scores_passed_through(graphml_files):
    mn_path, pr_path, base = graphml_files
    priority_df = pd.DataFrame({"feature_id": ["A", "B"], "triage_score": [0.8, 0.6]})
    rew_df = pd.DataFrame({"feature_id": ["A"], "rewiring_score": [0.75]})
    out = base / "out_scores"
    summary, _, __ = build_annotation_expansion(
        phenorewire_graphml=pr_path,
        specreboot_graphml=mn_path,
        priority_features_df=priority_df,
        rewiring_df=rew_df,
        hop_depth=1,
        min_cosine=0.70,
        output_dir=out,
    )
    row_a = summary[summary["feature_id"] == "A"].iloc[0]
    assert float(row_a["rewiring_score"]) == pytest.approx(0.75)
    assert float(row_a["triage_score"]) == pytest.approx(0.8)
