"""Connectivity-based selection: what it finds, what it must not find, and when it refuses."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from phenorewire.config import LeanConfig
from phenorewire.connectivity import connectivity_selection, differential_connectivity
from phenorewire.run import run

N_FEAT = 60
BLOCK = list(range(10))          # rewired: correlated in A, independent in B
SHIFT = list(range(10, 20))      # abundance shift only, correlation unchanged


def _block(n, k, rho, rng):
    return np.sqrt(rho) * rng.standard_normal((1, n)) + np.sqrt(1 - rho) * rng.standard_normal((k, n))


def _states(n, rng, *, rewire=True, shift=True):
    a = rng.standard_normal((N_FEAT, n))
    b = rng.standard_normal((N_FEAT, n))
    a[BLOCK] = _block(n, len(BLOCK), 0.85, rng)
    b[BLOCK] = _block(n, len(BLOCK), 0.0 if rewire else 0.85, rng)
    a[SHIFT] = _block(n, len(SHIFT), 0.6, rng)
    b[SHIFT] = _block(n, len(SHIFT), 0.6, rng) + (2.0 if shift else 0.0)
    for arr in (a, b):            # identical marginals for the rewired block
        arr[BLOCK] = (arr[BLOCK] - arr[BLOCK].mean(1, keepdims=True)) / arr[BLOCK].std(1, ddof=1, keepdims=True)
    return a + 10.0, b + 10.0


def _frame(m, prefix):
    return pd.DataFrame(m, index=[f"F{i:03d}" for i in range(m.shape[0])],
                        columns=[f"{prefix}{j:03d}" for j in range(m.shape[1])])


def test_score_is_zero_for_identical_states():
    rng = np.random.default_rng(0)
    a = rng.standard_normal((15, 12))
    assert np.allclose(differential_connectivity(a, a.copy(), 5), 0.0)


def test_selects_pure_rewiring_and_not_abundance_shift():
    rng = np.random.default_rng(1)
    a, b = _states(25, rng)
    table, summary = connectivity_selection(_frame(a, "A"), _frame(b, "B"), top_k=5,
                                            n_permutations=100, fdr_alpha=0.05, seed=0)
    selected = set(table.index[table["selected_connectivity"]])
    block = {f"F{i:03d}" for i in BLOCK}
    shifted = {f"F{i:03d}" for i in SHIFT}
    assert len(selected & block) >= 8
    assert not selected & shifted
    assert summary["n_features_selected"] == len(selected)


# ---------------------------------------------------------------------------
# Through the pipeline
# ---------------------------------------------------------------------------

def _write_dataset(tmp_path, n):
    rng = np.random.default_rng(7)
    a, b = _states(n, rng)
    X = np.hstack([a, b])
    samples = [f"S{j:03d}" for j in range(2 * n)]
    feat = pd.DataFrame({"feature_id": [f"F{i:03d}" for i in range(N_FEAT)],
                         "mz": 100 + np.arange(N_FEAT), "rt": 1 + 0.1 * np.arange(N_FEAT)})
    feat = pd.concat([feat, pd.DataFrame(np.exp(X), columns=[f"{s} Area" for s in samples])], axis=1)
    meta = pd.DataFrame({"sample": samples, "group": ["A"] * n + ["B"] * n})
    feat.to_csv(tmp_path / "f.csv", index=False)
    meta.to_csv(tmp_path / "m.csv", index=False)


def _config(tmp_path, tag, **extra):
    return LeanConfig(**{
        "DATA_MATRIX": str(tmp_path / "f.csv"), "METADATA": str(tmp_path / "m.csv"),
        "OUTDIR": str(tmp_path / f"out_{tag}"), "ANALYSIS_MODE": "phenotype",
        "GROUP_DEFINITION": {"A": {"group": "A"}, "B": {"group": "B"}},
        "PHENO_GROUP_REF": "A", "PHENO_GROUP_CASE": "B",
        "FEATURE_ID_COL": "feature_id", "MZ_COL": "mz", "RT_COL": "rt",
        "INTENSITY_REGEX": "Area$", "META_SAMPLE_COL": "sample",
        "N_PERMUTATIONS_PHENO": 200, "MIN_FEATURES_HARD_STOP": 3,
        "NETWORK_MIN_EDGES_HARD_STOP": 1, "NETWORK_MIN_NODES_HARD_STOP": 2,
        "CONNECTIVITY": {"top_k": 5, "n_permutations": 100},
        "NULL_MODEL": {"n_permutations": 50, "n_configuration_models": 5},
        **extra,
    })


def _connectivity_meta(outdir):
    report = json.loads((outdir / "report.json").read_text(encoding="utf-8"))
    return report["phenotype"]["selection"]["connectivity"]


def test_runs_automatically_and_reports_what_abundance_misses(tmp_path):
    # 25 per state: at 20 this 60-feature design recovers only part of the block
    # after BH (power, not a defect: see dev/benchmarks/null_benchmark.py).
    _write_dataset(tmp_path, 25)
    run(_config(tmp_path, "auto"))
    block = {f"F{i:03d}" for i in BLOCK}

    comp = tmp_path / "out_auto" / "A_vs_B"
    tested = pd.read_csv(comp / "details/selection/all_features_tested.csv", dtype={"feature_id": str})
    abundance_sel = tested[tested["significant"]]
    assert not set(abundance_sel["feature_id"]) & block, "the rewired block has no abundance signal"

    sel = pd.read_csv(comp / "selected_metabolites.csv", dtype={"feature_id": str})
    routes = dict(zip(sel["feature_id"], sel["selection_route"]))
    assert len([f for f in block if routes.get(f) == "connectivity"]) >= 6
    assert (comp / "details/selection/all_features_connectivity.csv").exists()

    meta = _connectivity_meta(tmp_path / "out_auto")
    assert meta["connectivity_run"] is True and meta["n_connectivity_only"] >= 6
    text = (tmp_path / "out_auto/REPORT.md").read_text(encoding="utf-8")
    assert "invisible to a differential-abundance analysis" in text


def test_abundance_only_below_operating_range_with_reason(tmp_path):
    _write_dataset(tmp_path, 8)
    run(_config(tmp_path, "small", N_PERMUTATIONS_PHENO=1000))
    meta = _connectivity_meta(tmp_path / "out_small")
    assert meta["connectivity_run"] is False and "at least 14 samples" in meta["reason"]
    sel = pd.read_csv(tmp_path / "out_small/A_vs_B/selected_metabolites.csv")
    assert "selection_route" not in sel.columns
    text = (tmp_path / "out_small/REPORT.md").read_text(encoding="utf-8")
    assert "Connectivity-based selection did not run" in text
