import numpy as np
import pandas as pd

from stats_utils import (
    adjust_pvalues_bh,
    legacy_fisher_z_transform,
    lookup_bh_adjusted,
    legacy_z_test,
    merge_filtered_unfiltered_chunk,
)
from tests.reference_legacy import (
    compare_filtered_unfiltered,
    exec_fish_z,
    fisher_z_transform,
    replace_with_mini_cor,
    z_test,
)


def test_legacy_fisher_and_z_match_original_functions():
    r = pd.Series([0.1, -0.5, 0.9999999999, np.nan])
    n = pd.Series([10, 20, 30, 12])
    expected_z = fisher_z_transform(r, n).to_numpy()
    actual_z = legacy_fisher_z_transform(r.to_numpy(), n.to_numpy())
    np.testing.assert_allclose(actual_z, expected_z, equal_nan=True)

    expected_p = z_test(expected_z, expected_z[::-1], n, n[::-1])
    actual_p = legacy_z_test(
        actual_z,
        actual_z[::-1],
        n.to_numpy(),
        n.to_numpy()[::-1],
    )
    np.testing.assert_allclose(actual_p, expected_p, equal_nan=True)


def test_chunk_merge_matches_original_dataframe_pipeline():
    filtered_corr = np.array([
        0.1,
        0.2000000004,
        np.nan,
        0.3,
        np.nan,
        0.8,
        -0.9,
        0.0,
    ])
    unfiltered_corr = np.array([
        0.1,
        0.2,
        0.4,
        np.nan,
        np.nan,
        0.2,
        -0.1,
        0.7,
    ])
    filtered_p = np.array([0.5, 0.4, np.nan, 0.03, np.nan, 0.001, 0.02, 0.9])
    unfiltered_p = np.array([0.5, 0.4, 0.02, np.nan, np.nan, 0.2, 0.8, 0.01])
    filtered_n = np.array([10, 10, 5, 6, 2, 20, 4, 30])
    unfiltered_n = np.array([10, 10, 5, 6, 2, 20, 4, 30])

    pairs = pd.DataFrame({
        "Substance 1": [f"A{i}" for i in range(len(filtered_corr))],
        "Substance 2": [f"B{i}" for i in range(len(filtered_corr))],
    })
    filtered_df = pairs.assign(
        Correlation=filtered_corr,
        **{"P-Value": filtered_p, "n": filtered_n},
    )
    unfiltered_df = pairs.assign(
        Correlation=unfiltered_corr,
        **{"P-Value": unfiltered_p, "n": unfiltered_n},
    )

    different = compare_filtered_unfiltered(filtered_df, unfiltered_df)
    significant = exec_fish_z(different)
    expected = replace_with_mini_cor(unfiltered_df, significant)

    actual_corr, actual_p = merge_filtered_unfiltered_chunk(
        filtered_corr,
        filtered_p,
        filtered_n,
        unfiltered_corr,
        unfiltered_p,
        unfiltered_n,
    )
    np.testing.assert_allclose(
        actual_corr,
        expected["Correlation"].to_numpy(),
        equal_nan=True,
    )
    np.testing.assert_allclose(
        actual_p,
        expected["P-Value"].to_numpy(),
        equal_nan=True,
    )


def _assert_low_memory_bh_matches_statsmodels(raw, tmp_path, alpha=0.05):
    from statsmodels.stats.multitest import multipletests

    metadata_path = tmp_path / "bh.json"
    result = adjust_pvalues_bh(
        raw,
        alpha=alpha,
        scan_chunk_size=3,
        metadata_path=metadata_path,
        show_progress=False,
    )

    valid = np.isfinite(raw)
    expected_full = np.full(raw.shape, np.nan, dtype=np.float64)
    if np.any(valid):
        expected_full[valid] = multipletests(
            raw[valid],
            alpha=alpha,
            method="fdr_bh",
        )[1]

    expected_rejected = valid & (expected_full <= alpha)
    actual_candidates = raw <= result["raw_cutoff"]
    np.testing.assert_array_equal(actual_candidates, expected_rejected)
    assert result["rejection_count"] == int(expected_rejected.sum())
    assert result["valid_count"] == int(valid.sum())

    actual_adjusted = lookup_bh_adjusted(
        raw[actual_candidates],
        result["sorted_significant_pvalues"],
        result["adjusted_significant"],
    )
    np.testing.assert_allclose(
        actual_adjusted,
        expected_full[actual_candidates],
        atol=1e-15,
        rtol=1e-13,
    )
    assert result["metadata"]["algorithm"] == (
        "exact_significant_prefix_fixed_point"
    )
    assert metadata_path.exists()


def test_low_memory_bh_matches_statsmodels_and_preserves_ties(tmp_path):
    raw = np.array(
        [0.01, np.nan, 0.2, 0.03, 0.001, 1.0, 0.03, 0.001, 0.05]
    )
    _assert_low_memory_bh_matches_statsmodels(raw, tmp_path)


def test_low_memory_bh_matches_statsmodels_randomized(tmp_path):
    rng = np.random.default_rng(1234)
    raw = rng.uniform(0.0, 1.0, size=2000)
    raw[:80] = rng.uniform(0.0, 1e-5, size=80)
    raw[100:110] = 0.0005
    raw[::113] = np.nan
    _assert_low_memory_bh_matches_statsmodels(raw, tmp_path)


def test_low_memory_bh_handles_no_valid_or_no_rejections(tmp_path):
    all_nan = np.full(12, np.nan)
    result = adjust_pvalues_bh(
        all_nan,
        scan_chunk_size=4,
        show_progress=False,
    )
    assert result["valid_count"] == 0
    assert result["rejection_count"] == 0
    assert result["sorted_significant_pvalues"].size == 0

    no_reject = np.linspace(0.2, 1.0, 50)
    result = adjust_pvalues_bh(
        no_reject,
        scan_chunk_size=7,
        show_progress=False,
    )
    assert result["valid_count"] == 50
    assert result["rejection_count"] == 0
