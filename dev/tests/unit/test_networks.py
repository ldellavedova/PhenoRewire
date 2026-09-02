import pandas as pd

from phenorewire.networks import (
    correlation_network,
    compute_network_metrics,
    compute_rewiring,
    summarize_rewiring,
    edge_id,
    check_louvain_stability,
)


def test_network_metrics_add_expected_columns() -> None:
    nodes = pd.DataFrame({"name": ["A", "B", "C"]}, index=["A", "B", "C"])
    edges = pd.DataFrame(
        {
            "node_a": ["A", "B"],
            "node_b": ["B", "C"],
            "r": [0.9, -0.8],
            "q": [0.01, 0.02],
            "sign": ["pos", "neg"],
        }
    )
    nodes_out, communities, summary = compute_network_metrics(nodes, edges)
    assert {"degree", "betweenness_centrality", "eigenvector_centrality", "community_id"}.issubset(nodes_out.columns)
    assert int(summary.iloc[0]["n_edges"]) == 2
    assert not communities.empty


def test_rewiring_summary_counts_edge_states() -> None:
    edges_a = pd.DataFrame(
        {
            "node_a": ["A", "A"],
            "node_b": ["B", "C"],
            "r": [0.7, 0.8],
            "q": [0.01, 0.01],
            "sign": ["pos", "pos"],
        }
    )
    edges_b = pd.DataFrame(
        {
            "node_a": ["A", "B"],
            "node_b": ["B", "C"],
            "r": [-0.7, 0.9],
            "q": [0.01, 0.01],
            "sign": ["neg", "pos"],
        }
    )
    _, _, _, union = compute_rewiring(edges_a, edges_b, label_a="ref", label_b="case")
    summary, nodes = summarize_rewiring(union, label_a="ref", label_b="case")
    assert int(summary.iloc[0]["n_shared_edges"]) == 1
    assert int(summary.iloc[0]["n_sign_switch_edges"]) == 1
    assert "rewiring_score" in nodes.columns


def test_correlation_network_relaxes_thresholds_when_base_is_too_strict() -> None:
    X = pd.DataFrame(
        [
            [1, 2, 3, 4, 5, 6, 7, 8],
            [1, 4, 2, 5, 3, 6, 8, 7],
            [8, 7, 6, 5, 4, 3, 2, 1],
            [8, 5, 7, 4, 6, 3, 2, 1],
        ],
        index=["A", "B", "C", "D"],
        columns=[f"s{i}" for i in range(8)],
        dtype=float,
    )
    _, edges, meta = correlation_network(
        X,
        abs_r_thr=0.9,
        fdr_alpha=0.01,
        adaptive_thresholds=True,
        min_edges=3,
        max_density=0.6,
        floor_abs_r_thr=0.5,
        ceil_fdr_alpha=0.2,
    )
    assert meta["adaptive_thresholds_used"] is True
    assert meta["applied_abs_r_thr"] < 0.9
    assert meta["selected_edge_count"] >= 3
    assert not edges.empty


def test_correlation_network_keeps_base_thresholds_when_signal_is_already_sufficient() -> None:
    X = pd.DataFrame(
        {
            "s0": [1, 1.2, 5, 5.1],
            "s1": [2, 2.1, 4, 4.2],
            "s2": [3, 3.2, 3, 3.1],
            "s3": [4, 4.1, 2, 2.2],
            "s4": [5, 5.0, 1, 1.1],
            "s5": [6, 6.1, 0, 0.2],
        },
        index=["A", "B", "C", "D"],
        dtype=float,
    )
    _, edges, meta = correlation_network(
        X,
        abs_r_thr=0.5,
        fdr_alpha=0.2,
        adaptive_thresholds=True,
        min_edges=1,
        max_density=1.0,
        floor_abs_r_thr=0.5,
        ceil_fdr_alpha=0.2,
    )
    assert meta["applied_abs_r_thr"] == 0.5
    assert meta["applied_fdr_alpha"] == 0.2
    assert meta["adaptive_thresholds_used"] is False
    assert meta["selected_edge_count"] == len(edges)


# -------------------------------------------------------------------
# edge_id
# -------------------------------------------------------------------

def test_edge_id_is_order_invariant() -> None:
    assert edge_id("A", "B") == edge_id("B", "A")
    assert edge_id("Z", "A") == edge_id("A", "Z")
    assert edge_id("X", "X") == "X__X"


# -------------------------------------------------------------------
# compute_rewiring corner cases
# -------------------------------------------------------------------

def _make_edges(pairs: list[tuple[str, str, float, str]]) -> pd.DataFrame:
    rows = [
        {"node_a": a, "node_b": b, "r": r, "q": 0.01, "sign": s}
        for a, b, r, s in pairs
    ]
    return pd.DataFrame(rows)


def test_rewiring_proportion_zero_when_identical_networks() -> None:
    edges = _make_edges([("A", "B", 0.8, "pos"), ("B", "C", 0.7, "pos")])
    _, _, _, union = compute_rewiring(edges, edges.copy(), label_a="ref", label_b="case")
    summary, _ = summarize_rewiring(union, label_a="ref", label_b="case")
    assert float(summary.iloc[0]["rewiring_proportion"]) == 0.0


def test_rewiring_proportion_one_when_no_shared_edges() -> None:
    edges_a = _make_edges([("A", "B", 0.8, "pos")])
    edges_b = _make_edges([("C", "D", 0.7, "pos")])
    _, _, _, union = compute_rewiring(edges_a, edges_b, label_a="ref", label_b="case")
    summary, _ = summarize_rewiring(union, label_a="ref", label_b="case")
    assert float(summary.iloc[0]["rewiring_proportion"]) == 1.0


def test_sign_switch_min_r_filters_marginal_correlations() -> None:
    """A sign switch between edges with |r| < min_r should NOT be flagged."""
    edges_a = _make_edges([("A", "B", 0.3, "pos"), ("C", "D", 0.8, "pos")])
    edges_b = _make_edges([("A", "B", -0.3, "neg"), ("C", "D", -0.8, "neg")])

    # Without filter → both flagged
    _, _, _, union_no_filter = compute_rewiring(
        edges_a, edges_b, label_a="ref", label_b="case", sign_switch_min_r=0.0
    )
    n_no_filter = int(union_no_filter["sign_switch"].fillna(False).sum())
    assert n_no_filter == 2

    # With min_r=0.5 → only the high-|r| pair (C-D) flagged
    _, _, _, union_filtered = compute_rewiring(
        edges_a, edges_b, label_a="ref", label_b="case", sign_switch_min_r=0.5
    )
    n_filtered = int(union_filtered["sign_switch"].fillna(False).sum())
    assert n_filtered == 1, (
        f"Expected 1 sign switch after min_r=0.5 filter, got {n_filtered}"
    )


def test_rewiring_manually_constructed_toy() -> None:
    """Verify set-level edge classification on a known toy example."""
    edges_a = _make_edges([
        ("A", "B", 0.9, "pos"),  # shared, no sign switch
        ("A", "C", 0.8, "pos"),  # REF-only
    ])
    edges_b = _make_edges([
        ("A", "B", -0.9, "neg"),  # shared, sign switch
        ("B", "C", 0.7, "pos"),   # CASE-only
    ])
    shared, only_a, only_b, union = compute_rewiring(
        edges_a, edges_b, label_a="ref", label_b="case", sign_switch_min_r=0.5
    )
    summary, nodes = summarize_rewiring(union, label_a="ref", label_b="case")

    assert int(summary.iloc[0]["n_shared_edges"]) == 1
    assert int(summary.iloc[0]["n_ref_only_edges"]) == 1
    assert int(summary.iloc[0]["n_case_only_edges"]) == 1
    assert int(summary.iloc[0]["n_sign_switch_edges"]) == 1
    assert "rewiring_score" in nodes.columns


# -------------------------------------------------------------------
# Louvain stability
# -------------------------------------------------------------------

def test_louvain_stability_perfect_for_trivial_graph() -> None:
    """A triangle always produces the same community — AMI should be 1.0."""
    edges = _make_edges([("A", "B", 0.9, "pos"), ("B", "C", 0.9, "pos"), ("A", "C", 0.9, "pos")])
    meta = check_louvain_stability(edges, n_runs=5, base_seed=42)
    assert meta["stable"] is True
    assert meta["ami_mean"] >= 0.99


def test_louvain_stability_empty_graph() -> None:
    meta = check_louvain_stability(None, n_runs=5, base_seed=42)
    assert meta["stable"] is True
    assert meta["warning"] is None
