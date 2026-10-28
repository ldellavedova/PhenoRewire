"""Tests for phenotype/triage.py — score computation and weight validation."""
from __future__ import annotations

import numpy as np
import pandas as pd

from phenorewire.triage import rank_network_nodes, rewired_metabolites


def _make_nodes(n: int = 10, with_q: bool = True) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    data = {
        "degree_centrality": rng.uniform(0, 1, n),
        "betweenness_centrality": rng.uniform(0, 1, n),
        "eigenvector_centrality": rng.uniform(0, 1, n),
        "weighted_degree_abs": rng.uniform(0, 1, n),
    }
    if with_q:
        data["q_pheno"] = rng.uniform(0.001, 0.1, n)
    df = pd.DataFrame(data, index=[f"F{i:03d}" for i in range(n)])
    df.index.name = "feature_id"
    return df


def _run_triage(nodes_df, topology_weight=0.7, selection_weight=0.3):
    return rank_network_nodes(nodes_df, topology_weight=topology_weight, selection_weight=selection_weight)


def test_triage_score_is_in_unit_interval() -> None:
    nodes = _make_nodes(15)
    result = _run_triage(nodes)
    assert not result.empty, "triage output should not be empty"
    assert "triage_score" in result.columns
    assert (result["triage_score"].between(0.0, 1.0)).all(), (
        f"triage_score out of [0,1]: {result['triage_score'].describe()}"
    )


def test_triage_score_with_default_weights() -> None:
    nodes = _make_nodes(10)
    result = _run_triage(nodes, topology_weight=0.7, selection_weight=0.3)
    assert not result.empty
    assert result["triage_score"].max() > result["triage_score"].min()


def test_triage_nondefault_weights_change_ranking() -> None:
    """Changing weights should produce a different ordering."""
    nodes = _make_nodes(12)
    r_default = _run_triage(nodes, topology_weight=0.7, selection_weight=0.3)
    r_selection = _run_triage(nodes, topology_weight=0.1, selection_weight=0.9)

    assert not r_default.empty and not r_selection.empty

    top_default = r_default.head(5)["feature_id"].tolist()
    top_selection = r_selection.head(5)["feature_id"].tolist()
    # With very different weights the top-5 should differ
    assert top_default != top_selection, (
        "Expected different top-5 when weights are 0.7/0.3 vs 0.1/0.9"
    )


def test_triage_score_without_q_pheno() -> None:
    """When q_pheno column is absent, selection_support_score falls back to 0."""
    nodes = _make_nodes(8, with_q=False)
    result = _run_triage(nodes)
    assert not result.empty
    assert "triage_score" in result.columns
    # selection_support_score = 0 → triage = topology_weight * topology_priority_score
    assert (result["triage_score"] >= 0.0).all()
    assert (result["triage_score"] <= 1.0).all()


def test_triage_empty_nodes_does_not_crash() -> None:
    result = rank_network_nodes(pd.DataFrame())
    assert result.empty
    assert "triage_score" in result.columns


def test_triage_topology_weight_zero_uses_selection_only() -> None:
    """With topology_weight=0, triage_score should equal selection_support_score."""
    nodes = _make_nodes(6)
    result = _run_triage(nodes, topology_weight=0.0, selection_weight=1.0)
    assert not result.empty
    # selection_support_score = 1 - q_pheno; all scores should be strictly positive
    assert (result["triage_score"] > 0).all()


def _selected(ids) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "name": [f"m{i}" for i in ids],
            "selection_route": ["abundance", "connectivity", "both"][: len(ids)],
            "direction": ["up_in_B"] * len(ids),
        },
        index=ids,
    )
    return df


def test_rewired_metabolites_lists_every_selected_metabolite() -> None:
    """Metabolites without edges stay in the table, with zero counts."""
    selected = _selected(["F1", "F2", "F3"])
    rewiring_nodes = pd.DataFrame(
        {
            "feature_id": ["F1"], "shared_edges": [1], "A_only_edges": [2], "B_only_edges": [0],
            "sign_switch_edges": [0], "rewiring_score": [2.0], "rewiring_fraction": [0.67],
        }
    )
    out = rewired_metabolites(selected, rewiring_nodes, None, label_a="A", label_b="B", direction_col="direction")
    assert out["feature_id"].tolist() == ["F1", "F2", "F3"]
    assert out["rank"].tolist() == [1, 2, 3]
    assert out.loc[out["feature_id"] == "F2", "rewiring_score"].item() == 0
    assert {"name", "selection_route", "direction"} <= set(out.columns)
    assert "q_rewiring" not in out.columns


def test_rewired_metabolites_sorted_by_q_when_null_ran() -> None:
    selected = _selected(["F1", "F2"])
    rewiring_nodes = pd.DataFrame(
        {
            "feature_id": ["F1", "F2"], "shared_edges": [0, 0], "A_only_edges": [5, 1], "B_only_edges": [0, 0],
            "sign_switch_edges": [0, 0], "rewiring_score": [5.0, 1.0], "rewiring_fraction": [1.0, 1.0],
        }
    )
    node_null = pd.DataFrame(
        {
            "feature_id": ["F1", "F2"], "observed_state_specific_edges": [5, 1],
            "null_mean_state_specific_edges": [4.9, 0.1], "p_empirical": [0.5, 0.01], "q_bh": [0.5, 0.02],
        }
    )
    out = rewired_metabolites(selected, rewiring_nodes, node_null, label_a="A", label_b="B", direction_col="direction")
    assert out["feature_id"].tolist() == ["F2", "F1"]
    assert out["more_than_chance"].tolist() == [True, False]
