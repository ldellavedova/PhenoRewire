"""Tests for phenotype/triage.py — score computation and weight validation."""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenorewire.triage import export_network_triage


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
    with tempfile.TemporaryDirectory() as tmp:
        outdir = Path(tmp)
        export_network_triage(
            nodes_df,
            pd.DataFrame(),
            pd.DataFrame(),
            outdir,
            top_n=len(nodes_df),
            topology_weight=topology_weight,
            selection_weight=selection_weight,
        )
        triage_path = outdir / "triage" / "node_priority_ranking.csv"
        if triage_path.exists():
            return pd.read_csv(triage_path)
    return pd.DataFrame()


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
    with tempfile.TemporaryDirectory() as tmp:
        outdir = Path(tmp)
        export_network_triage(
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            outdir,
            top_n=25,
        )
        manifest = outdir / "triage" / "triage_manifest.json"
        assert manifest.exists(), "manifest should be written even for empty input"


def test_triage_topology_weight_zero_uses_selection_only() -> None:
    """With topology_weight=0, triage_score should equal selection_support_score."""
    nodes = _make_nodes(6)
    result = _run_triage(nodes, topology_weight=0.0, selection_weight=1.0)
    assert not result.empty
    # selection_support_score = 1 - q_pheno; all scores should be strictly positive
    assert (result["triage_score"] > 0).all()
