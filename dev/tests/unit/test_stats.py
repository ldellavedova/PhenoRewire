import numpy as np
import pytest
from scipy.stats import spearmanr

from phenorewire.stats import adaptive_fdr_selection, spearman_permutation_test
from phenorewire.phenotype_selection import phenotype_selection


def test_spearman_matches_scipy_on_observed_correlations() -> None:
    mat = np.array(
        [
            [1, 2, 3, 4, 5, 6],
            [6, 5, 4, 3, 2, 1],
            [1, 1, 1, 1, 1, 1],
        ],
        dtype=float,
    )
    target = np.array([0, 0, 0, 1, 1, 1], dtype=float)
    r_obs, p_emp = spearman_permutation_test(mat, target, n_permutations=20, seed=42)
    assert np.isclose(r_obs[0], spearmanr(mat[0], target).statistic)
    assert np.isclose(r_obs[1], spearmanr(mat[1], target).statistic)
    assert np.isnan(r_obs[2])
    assert np.all((p_emp >= 0.0) & (p_emp <= 1.0))


def test_same_seed_reproduces_identical_results() -> None:
    mat = np.array(
        [
            [1, 3, 2, 4, 5, 8],
            [8, 5, 4, 3, 2, 1],
        ],
        dtype=float,
    )
    target = np.array([0, 1, 0, 1, 0, 1], dtype=float)
    r_a, p_a = spearman_permutation_test(mat, target, n_permutations=25, seed=7)
    r_b, p_b = spearman_permutation_test(mat, target, n_permutations=25, seed=7)
    assert np.array_equal(r_a, r_b, equal_nan=True)
    assert np.array_equal(p_a, p_b, equal_nan=True)


def test_rows_with_nan_fall_back_and_do_not_poison_other_features() -> None:
    """A NaN in one feature must not change any other feature's result."""
    rng = np.random.default_rng(3)
    mat = rng.normal(size=(6, 12))
    target = np.repeat([0.0, 1.0], 6)
    mat[0] += target * 2.0

    r_clean, p_clean = spearman_permutation_test(mat, target, n_permutations=50, seed=5)

    holed = mat.copy()
    holed[3, 2] = np.nan
    r_holed, p_holed = spearman_permutation_test(holed, target, n_permutations=50, seed=5)

    untouched = [0, 1, 2, 4, 5]
    assert np.array_equal(r_clean[untouched], r_holed[untouched], equal_nan=True)
    assert np.array_equal(p_clean[untouched], p_holed[untouched], equal_nan=True)
    assert np.isfinite(r_holed[3])


def test_constant_feature_returns_nan_and_p_one() -> None:
    mat = np.array([[1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
    target = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    r_obs, p_emp = spearman_permutation_test(mat, target, n_permutations=20, seed=1)
    assert np.isnan(r_obs[0])
    assert p_emp[0] == 1.0
    assert np.isfinite(r_obs[1])


def test_target_with_nan_is_rejected() -> None:
    mat = np.ones((2, 6))
    mat[1] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    target = np.array([0.0, 1.0, np.nan, 1.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="NaN"):
        spearman_permutation_test(mat, target, n_permutations=20, seed=1)


def test_adaptive_fdr_selection_relaxes_alpha_when_signal_is_borderline() -> None:
    qvals = np.array([0.0526, 0.0526, 0.0526, 0.08, 0.2])
    selected, meta = adaptive_fdr_selection(
        qvals,
        base_alpha=0.01,
        adaptive=True,
        min_selected=3,
        alpha_ceiling=0.1,
    )
    assert meta["adaptive_fdr_used"] is True
    assert np.isclose(meta["applied_fdr_alpha"], 0.1)
    assert int(selected.sum()) >= 3


def test_adaptive_fdr_selection_keeps_base_alpha_when_already_sufficient() -> None:
    qvals = np.array([0.001, 0.01, 0.02, 0.2])
    selected, meta = adaptive_fdr_selection(
        qvals,
        base_alpha=0.05,
        adaptive=True,
        min_selected=2,
        alpha_ceiling=0.1,
    )
    assert meta["adaptive_fdr_used"] is False
    assert np.isclose(meta["applied_fdr_alpha"], 0.05)
    assert int(selected.sum()) == 3


def test_empirical_p_low_for_strong_correlation() -> None:
    """A perfectly monotone feature should yield empirical p < 0.05."""
    rng = np.random.default_rng(1)
    n = 20
    y = np.array([0] * 10 + [1] * 10, dtype=float)
    # Feature that perfectly tracks y
    x_signal = y + rng.normal(0, 0.01, n)
    mat = x_signal.reshape(1, -1)
    r_obs, p_emp = spearman_permutation_test(mat, y, n_permutations=500, seed=99)
    assert p_emp[0] < 0.05, f"Strong correlation should have p<0.05, got {p_emp[0]}"


def test_empirical_p_not_low_for_noise() -> None:
    """Pure noise features should not systematically yield p < 0.05."""
    rng = np.random.default_rng(2)
    n_feat, n_samp = 50, 20
    y = np.array([0] * 10 + [1] * 10, dtype=float)
    mat = rng.standard_normal((n_feat, n_samp))
    r_obs, p_emp = spearman_permutation_test(mat, y, n_permutations=200, seed=7)
    # Expect very few false positives at p<0.05 under the null
    n_significant = int((p_emp < 0.05).sum())
    assert n_significant < 10, f"Too many false positives in noise: {n_significant}/50"


def test_phenotype_selection_hard_stop_raises_on_no_signal(tmp_path) -> None:
    """Pipeline should raise when selected features < min_features_hard_stop."""
    rng = np.random.default_rng(3)
    n_feat, n_samp = 10, 20
    y = np.array([0] * 10 + [1] * 10, dtype=int)
    # Pure noise: nothing should be selected
    mat = rng.standard_normal((n_feat, n_samp))
    X = __import__("pandas").DataFrame(
        mat,
        index=[f"F{i}" for i in range(n_feat)],
        columns=[f"S{i}" for i in range(n_samp)],
    )
    import pandas as pd
    feat_anno = pd.DataFrame(index=X.index)
    feat_anno.index.name = "feature_id"

    with pytest.raises(ValueError, match="hard stop"):
        phenotype_selection(
            X, y, feat_anno, tmp_path,
            n_permutations=50,
            fdr_alpha=0.001,         # extremely strict → 0 selected
            adaptive_fdr=False,
            min_selected_features=5,
            fdr_alpha_ceiling=0.001,
            seed=42,
            min_features_hard_stop=5,
        )


def test_adaptive_fdr_off_by_default_does_not_relax() -> None:
    """With adaptive=False, applied alpha must equal base alpha."""
    qvals = np.array([0.15, 0.2, 0.25, 0.3])
    selected, meta = adaptive_fdr_selection(
        qvals,
        base_alpha=0.05,
        adaptive=False,
        min_selected=3,
        alpha_ceiling=0.5,
    )
    assert meta["adaptive_fdr_used"] is False
    assert np.isclose(meta["applied_fdr_alpha"], 0.05)
    assert int(selected.sum()) == 0


# ---------------------------------------------------------------------------
# Permutation streams
# ---------------------------------------------------------------------------

def _legacy_shared_stream(X, y, n_perm, seed):
    """The pre-release implementation, verbatim in behaviour: one generator,
    one scipy.spearmanr call per permutation, constant rows skipped."""
    from scipy.stats import spearmanr

    rng = np.random.default_rng(seed)
    r_obs = np.full(X.shape[0], np.nan)
    p = np.ones(X.shape[0])
    for i in range(X.shape[0]):
        x = X[i]
        if np.nanstd(x) == 0.0:
            continue
        r, _ = spearmanr(x, y)
        r_obs[i] = r
        rp = np.array([abs(spearmanr(x, rng.permutation(y))[0]) for _ in range(n_perm)])
        p[i] = (np.sum(rp >= abs(r)) + 1) / (n_perm + 1)
    return r_obs, p


def test_shared_stream_reproduces_pre_release_scheme():
    rng = np.random.default_rng(11)
    X = rng.standard_normal((12, 10))
    X[4] = 3.0                                   # constant row draws nothing
    y = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=float)
    r_new, p_new = spearman_permutation_test(X, y, n_permutations=50, seed=42, stream="shared")
    r_old, p_old = _legacy_shared_stream(X, y, 50, 42)
    np.testing.assert_allclose(r_new, r_old, atol=1e-12, equal_nan=True)
    np.testing.assert_array_equal(p_new, p_old)


def test_streams_agree_on_r_and_differ_only_in_draws():
    rng = np.random.default_rng(12)
    X = rng.standard_normal((30, 16))
    y = np.repeat([0.0, 1.0], 8)
    r1, p1 = spearman_permutation_test(X, y, n_permutations=200, seed=1, stream="per_feature")
    r2, p2 = spearman_permutation_test(X, y, n_permutations=200, seed=1, stream="shared")
    np.testing.assert_allclose(r1, r2, atol=1e-12)
    assert np.max(np.abs(p1 - p2)) < 0.15


def test_unknown_stream_rejected():
    with pytest.raises(ValueError, match="stream"):
        spearman_permutation_test(np.ones((2, 5)), np.arange(5.0), n_permutations=10, seed=0, stream="bogus")
