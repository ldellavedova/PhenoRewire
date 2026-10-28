"""Null models: identical filtering to the real networks, and calibration."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenorewire.config import NullModelConfig
from phenorewire.networks import correlation_network
from phenorewire.nulls import (
    edge_index_set,
    modularity_null,
    network_edge_null,
    rewiring_label_null,
)
from phenorewire.run import operating_range


def _block(n: int, k: int, rho: float, rng) -> np.ndarray:
    """k features x n samples; equicorrelated at rho."""
    shared = rng.standard_normal((1, n))
    return np.sqrt(rho) * shared + np.sqrt(1 - rho) * rng.standard_normal((k, n))


def _frame(mat: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        mat,
        index=[f"F{i:03d}" for i in range(mat.shape[0])],
        columns=[f"S{j:03d}" for j in range(mat.shape[1])],
    )


# ---------------------------------------------------------------------------
# Same filter as the observed networks
# ---------------------------------------------------------------------------

def test_edge_index_set_matches_correlation_network():
    rng = np.random.default_rng(1)
    mat = np.vstack([_block(20, 8, 0.9, rng), rng.standard_normal((10, 20))])
    mat[3] = 5.0                     # constant row, dropped by both
    mat[5, 2] = np.nan               # one missing value, removes the feature in both
    X = _frame(mat)

    _, edges_df, _ = correlation_network(X, abs_r_thr=0.6, fdr_alpha=0.05, adaptive_thresholds=False)
    expected = {tuple(sorted((a, b))) for a, b in zip(edges_df["node_a"], edges_df["node_b"])}

    feats = X.index.tolist()
    got = {tuple(sorted((feats[i], feats[j]))) for i, j in edge_index_set(mat, 0.6, 0.05)}
    assert expected, "fixture should produce edges"
    assert got == expected


# ---------------------------------------------------------------------------
# Edge null
# ---------------------------------------------------------------------------

def test_edge_null_flags_structure_and_not_noise():
    rng = np.random.default_rng(2)
    structured = _frame(np.vstack([_block(20, 10, 0.8, rng), rng.standard_normal((10, 20))]))
    noise = _frame(rng.standard_normal((20, 20)))

    s = network_edge_null(structured, abs_r_thr=0.6, fdr_alpha=0.05, n_permutations=50, seed=0)
    n = network_edge_null(noise, abs_r_thr=0.6, fdr_alpha=0.05, n_permutations=50, seed=0)

    assert s["observed_edges"] > 20
    assert s["p_empirical"] <= 1 / 51 + 1e-12
    assert n["p_empirical"] > 0.05


def test_edge_null_is_reproducible():
    rng = np.random.default_rng(3)
    X = _frame(np.vstack([_block(15, 6, 0.7, rng), rng.standard_normal((6, 15))]))
    a = network_edge_null(X, abs_r_thr=0.6, fdr_alpha=0.05, n_permutations=30, seed=7)
    b = network_edge_null(X, abs_r_thr=0.6, fdr_alpha=0.05, n_permutations=30, seed=7)
    assert a == b


def test_modularity_null_two_cliques_beats_configuration_model():
    rows = []
    for block in (range(0, 6), range(6, 12)):
        for i in block:
            for j in block:
                if i < j:
                    rows.append((f"F{i:02d}", f"F{j:02d}", 0.9, 0.001, "pos"))
    rows.append(("F00", "F06", 0.9, 0.001, "pos"))
    edges = pd.DataFrame(rows, columns=["node_a", "node_b", "r", "q", "sign"])
    out = modularity_null(edges, n_reps=20, seed=0)
    assert out["modularity_unweighted"] > 0.4
    assert out["delta_modularity"] > 0.2


# ---------------------------------------------------------------------------
# Rewiring null
# ---------------------------------------------------------------------------

def _two_states(n: int, rng, *, shift: float, rewire: bool):
    """30 features; block F000-F009 correlated in A, and in B unless ``rewire``."""
    a = np.vstack([_block(n, 10, 0.8, rng), rng.standard_normal((20, n))])
    b_block = rng.standard_normal((10, n)) if rewire else _block(n, 10, 0.8, rng)
    b = np.vstack([b_block, rng.standard_normal((20, n))])
    b[:10] += shift
    return _frame(a), _frame(b).set_axis([f"T{j:03d}" for j in range(n)], axis=1)


def test_rewiring_null_detects_rewiring():
    rng = np.random.default_rng(4)
    A, B = _two_states(20, rng, shift=1.0, rewire=True)
    summary, nodes = rewiring_label_null(
        A, B, thresholds_a=(0.6, 0.05), thresholds_b=(0.6, 0.05), n_permutations=100, seed=0
    )
    assert summary["p_empirical"] <= 0.02
    assert summary["obs_over_null"] > 2
    top = set(nodes.loc[nodes["q_bh"] <= 0.05, "feature_id"])
    assert top and top <= {f"F{i:03d}" for i in range(10)}


def test_rewiring_null_is_calibrated_when_only_abundance_changes():
    """The case raw state-specific edge counts get wrong: a block shifts its mean
    but keeps its correlation.  Centring within state must keep the null valid,
    so p-values stay roughly uniform rather than piling up near zero."""
    rng = np.random.default_rng(5)
    pvals = []
    for _ in range(30):
        A, B = _two_states(16, rng, shift=1.5, rewire=False)
        summary, _ = rewiring_label_null(
            A, B, thresholds_a=(0.6, 0.05), thresholds_b=(0.6, 0.05), n_permutations=30, seed=1
        )
        pvals.append(summary["p_empirical"])
    pvals = np.array(pvals)
    assert np.mean(pvals <= 0.05) <= 0.15
    assert 0.3 <= np.median(pvals) <= 0.75


def test_rewiring_null_requires_matching_features():
    rng = np.random.default_rng(6)
    A, B = _two_states(12, rng, shift=0.0, rewire=False)
    with pytest.raises(ValueError, match="same features"):
        rewiring_label_null(
            A, B.iloc[::-1], thresholds_a=(0.6, 0.05), thresholds_b=(0.6, 0.05), n_permutations=20, seed=0
        )


# ---------------------------------------------------------------------------
# Config and operating range
# ---------------------------------------------------------------------------

def test_null_model_config_defaults_and_validation():
    cfg = NullModelConfig()
    assert cfg.enabled is True and cfg.n_permutations == 200
    with pytest.raises(ValueError):
        NullModelConfig(n_permutations=10)
    with pytest.raises(ValueError):
        NullModelConfig(n_configuration_models=-1)


@pytest.mark.parametrize(
    "n_a, n_b, status",
    [(14, 30, "within"), (10, 13, "reduced"), (5, 20, "below"), (9, 9, "below")],
)
def test_operating_range_uses_the_smaller_state(n_a, n_b, status):
    assert operating_range(n_a, n_b)["status"] == status
