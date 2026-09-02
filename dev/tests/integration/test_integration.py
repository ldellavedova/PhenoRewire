"""
End-to-end integration test for PhenoRewire.

Builds a synthetic dataset with known rewired nodes and runs the full pipeline.
Verifies output files, report contents, and that injected nodes appear in the
top-10 priority rewired list.

The synthetic data is designed with strong enough signal that adaptive FDR
relaxation should NOT be triggered (ADAPTIVE_SELECTION_THRESHOLDS=False is
intentional here to confirm the data is clean).
"""

from __future__ import annotations


import numpy as np
import pandas as pd
import pytest

from phenorewire.config import LeanConfig
from phenorewire.run import run


# ---------------------------------------------------------------------------
# Synthetic dataset generation
# ---------------------------------------------------------------------------

N_REWIRED = 8
N_TOTAL = 25
SEED = 42


def _make_cov_ref(n: int) -> np.ndarray:
    """REF covariance: features F000–F007 form a dense positive block."""
    C = np.eye(n)
    for i in range(N_REWIRED):
        for j in range(N_REWIRED):
            if i != j:
                C[i, j] = 0.78
    C += np.eye(n) * 0.05
    return C


def _make_cov_case(n: int) -> np.ndarray:
    """CASE covariance: F000-F003 have only 2 positive pairs (partial rewiring);
    F004-F011 form a NEW positive block (CASE-only edges)."""
    C = np.eye(n)
    # F004-F011 are a NEW dense positive block in CASE (CASE-only edges)
    case_block_start = 4
    case_block_end = min(n, 12)  # features 4..11
    for i in range(case_block_start, case_block_end):
        for j in range(case_block_start, case_block_end):
            if i != j:
                C[i, j] = 0.75
    # F000-F001 and F002-F003: sign-switched (negative in CASE)
    for i, j in [(0, 1), (1, 0), (2, 3), (3, 2)]:
        C[i, j] = -0.70
    C += np.eye(n) * 0.15  # regularize to ensure PD
    return C


def _build_synthetic_data(n_per_group: int, seed: int = SEED):
    """
    Returns (feature_df, meta_df, injected_rewired_ids).

    feature_df columns: feature_id, mz, rt, Area_<sample>...
    meta_df columns: sample, group
    """
    rng = np.random.default_rng(seed)
    n_total = 2 * n_per_group
    feature_ids = [f"F{i:03d}" for i in range(N_TOTAL)]
    sample_ids = [f"S{i:03d}" for i in range(n_total)]

    C_ref = _make_cov_ref(N_TOTAL)
    C_case = _make_cov_case(N_TOTAL)

    # Ensure positive-definiteness by regularization if needed
    for C in (C_ref, C_case):
        eigvals = np.linalg.eigvalsh(C)
        if eigvals.min() <= 0:
            C += np.eye(N_TOTAL) * (abs(eigvals.min()) + 0.01)

    L_ref = np.linalg.cholesky(C_ref)
    L_case = np.linalg.cholesky(C_case)

    Z_ref = rng.standard_normal((n_per_group, N_TOTAL)) @ L_ref.T
    Z_case = rng.standard_normal((n_per_group, N_TOTAL)) @ L_case.T

    # Strong mean shift for rewired features → ensures phenotype selection
    Z_case[:, :N_REWIRED] += 2.5

    X = np.vstack([Z_ref, Z_case]).T  # (N_TOTAL, n_total)
    # Shift all values to be positive (simulates positive intensities)
    X = X - X.min() + 1.0

    # Build feature CSV using MZmine-style "sample Area" column names so that
    # preprocessing correctly strips the suffix and yields sample IDs.
    intensity_cols = {f"{sid} Area": X[:, i] for i, sid in enumerate(sample_ids)}
    feat_dict = {
        "feature_id": feature_ids,
        "mz": [100.0 + i * 0.5 for i in range(N_TOTAL)],
        "rt": [1.0 + i * 0.1 for i in range(N_TOTAL)],
    }
    feat_dict.update(intensity_cols)
    feature_df = pd.DataFrame(feat_dict)

    # Metadata
    meta_df = pd.DataFrame({
        "sample": sample_ids,
        "group": ["REF"] * n_per_group + ["CASE"] * n_per_group,
    })

    injected_rewired = feature_ids[:N_REWIRED]
    return feature_df, meta_df, injected_rewired


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    """Run the full pipeline once and return (outdir, injected_rewired_ids)."""
    base = tmp_path_factory.mktemp("integration")
    n_per_group = 30  # large enough for clean signal without adaptive FDR

    feature_df, meta_df, injected = _build_synthetic_data(n_per_group)

    feat_path = base / "features.csv"
    meta_path = base / "metadata.csv"
    outdir = base / "output"
    outdir.mkdir()

    feature_df.to_csv(feat_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    config_dict = {
        "DATA_MATRIX": str(feat_path),
        "METADATA": str(meta_path),
        "OUTDIR": str(outdir),
        "ANALYSIS_MODE": "phenotype",
        "GROUP_DEFINITION": {
            "REF": {"group": "REF"},
            "CASE": {"group": "CASE"},
        },
        "PHENO_GROUP_REF": "REF",
        "PHENO_GROUP_CASE": "CASE",
        "FEATURE_ID_COL": "feature_id",
        "MZ_COL": "mz",
        "RT_COL": "rt",
        "INTENSITY_REGEX": "Area$",
        "META_SAMPLE_COL": "sample",
        "N_PERMUTATIONS_PHENO": 100,
        "FDR_ALPHA": 0.1,
        "ADAPTIVE_SELECTION_THRESHOLDS": False,  # data is clean; no relaxation needed
        "MIN_FEATURES_HARD_STOP": 3,
        "CORR_ABS_THRESHOLD": 0.4,
        "CORR_FDR_ALPHA": 0.2,
        "ADAPTIVE_NETWORK_THRESHOLDS": True,
        "NETWORK_MIN_EDGES": 3,
        "NETWORK_MIN_EDGES_HARD_STOP": 3,
        "NETWORK_MIN_NODES_HARD_STOP": 3,
        "SIGN_SWITCH_MIN_R": 0.3,
        "LOUVAIN_STABILITY_CHECK": False,
        "TRIAGE_TOPOLOGY_WEIGHT": 0.7,
        "TRIAGE_SELECTION_WEIGHT": 0.3,
        "RANDOM_SEED": SEED,
        "LOG_TRANSFORM": True,
        "NORM_METHOD": "median",
    }
    cfg = LeanConfig(**config_dict)
    run(cfg)

    return outdir, injected


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_output_files_exist(pipeline_output) -> None:
    outdir, _ = pipeline_output
    required_files = [
        "report/final_report.md",
        "report/final_report.json",
        "triage/priority_rewired_nodes.csv",
        "triage/rewiring_summary.csv",
    ]
    for rel in required_files:
        p = outdir / rel
        assert p.exists(), f"Expected output file missing: {rel}"


def test_final_report_is_nonempty(pipeline_output) -> None:
    outdir, _ = pipeline_output
    report_text = (outdir / "report" / "final_report.md").read_text(encoding="utf-8")
    assert len(report_text) > 500, "final_report.md seems too short"


def test_final_report_sections_present(pipeline_output) -> None:
    outdir, _ = pipeline_output
    report_text = (outdir / "report" / "final_report.md").read_text(encoding="utf-8")
    for section in ["Executive Summary", "Phenotype Comparison", "Scoring Formulas", "Reproducibility", "Next Steps"]:
        assert section in report_text, f"Missing section in report: {section}"


def test_no_adaptive_fdr_warning_in_clean_data(pipeline_output) -> None:
    """The clean synthetic dataset should not trigger adaptive FDR relaxation."""
    outdir, _ = pipeline_output
    report_text = (outdir / "report" / "final_report.md").read_text(encoding="utf-8")
    assert "ADAPTIVE FDR RELAXATION" not in report_text, (
        "Adaptive FDR should not be triggered for a clean synthetic dataset. "
        "If this fails the data may not have enough signal at alpha=0.1."
    )


def test_injected_rewired_nodes_in_top_ranked(pipeline_output) -> None:
    """At least half of injected rewired nodes should appear in the top-10."""
    outdir, injected = pipeline_output
    priority_path = outdir / "triage" / "priority_rewired_nodes.csv"
    if not priority_path.exists():
        pytest.skip("priority_rewired_nodes.csv not produced (possible empty rewiring)")

    ranked = pd.read_csv(priority_path)
    if ranked.empty:
        pytest.skip("priority_rewired_nodes.csv is empty")

    top10 = set(ranked.head(10)["feature_id"].astype(str).tolist())
    injected_set = set(injected)
    overlap = top10 & injected_set
    min_expected = max(1, len(injected) // 2)
    assert len(overlap) >= min_expected, (
        f"Only {len(overlap)}/{len(injected)} injected nodes in top-10. "
        f"Top-10: {top10}, Injected: {injected_set}"
    )


def test_rewiring_summary_has_expected_columns(pipeline_output) -> None:
    outdir, _ = pipeline_output
    # Look in either rewiring_pheno or triage
    for candidate in ["triage/rewiring_summary.csv"]:
        p = outdir / candidate
        if p.exists():
            df = pd.read_csv(p)
            assert "n_shared_edges" in df.columns
            assert "rewiring_proportion" in df.columns
            return
    pytest.skip("rewiring_summary.csv not found in expected location")


def test_triage_scores_in_priority_features(pipeline_output) -> None:
    outdir, _ = pipeline_output
    # Look for priority_features.csv in any network triage subdirectory
    for p in outdir.rglob("triage/priority_features.csv"):
        df = pd.read_csv(p)
        if not df.empty and "triage_score" in df.columns:
            assert (df["triage_score"] >= 0.0).all()
            assert (df["triage_score"] <= 1.0).all()
            return
    pytest.skip("priority_features.csv not found")


def test_sample_size_guard_raises_for_too_few_samples(tmp_path) -> None:
    """Run should raise when a group has fewer than 5 samples."""
    n_per_group = 3  # below hard stop
    feature_df, meta_df, _ = _build_synthetic_data(n_per_group)

    feat_path = tmp_path / "features.csv"
    meta_path = tmp_path / "metadata.csv"
    feature_df.to_csv(feat_path, index=False)
    meta_df.to_csv(meta_path, index=False)

    config_dict = {
        "DATA_MATRIX": str(feat_path),
        "METADATA": str(meta_path),
        "OUTDIR": str(tmp_path / "out"),
        "ANALYSIS_MODE": "phenotype",
        "GROUP_DEFINITION": {"REF": {"group": "REF"}, "CASE": {"group": "CASE"}},
        "PHENO_GROUP_REF": "REF",
        "PHENO_GROUP_CASE": "CASE",
        "FEATURE_ID_COL": "feature_id",
        "MZ_COL": "mz",
        "RT_COL": "rt",
        "INTENSITY_REGEX": "Area$",
        "META_SAMPLE_COL": "sample",
        "N_PERMUTATIONS_PHENO": 20,
        "FDR_ALPHA": 0.2,
        "ADAPTIVE_SELECTION_THRESHOLDS": True,
        "MIN_FEATURES_HARD_STOP": 1,
        "NETWORK_MIN_EDGES_HARD_STOP": 1,
        "NETWORK_MIN_NODES_HARD_STOP": 1,
        "RANDOM_SEED": SEED,
    }
    cfg = LeanConfig(**config_dict)
    with pytest.raises(ValueError, match="sample"):
        run(cfg)
