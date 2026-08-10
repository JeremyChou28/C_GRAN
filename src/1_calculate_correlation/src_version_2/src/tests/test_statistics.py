import numpy as np
import pandas as pd

from stats_utils import (
    adjust_pvalues_bh,
    legacy_fisher_z_transform,
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


def test_bh_matches_statsmodels_and_preserves_nan(tmp_path):
    raw = np.array([0.01, np.nan, 0.2, 0.03, 0.001, 1.0, 0.03])
    result = adjust_pvalues_bh(raw, tmp_path / "adjusted.dat")

    from statsmodels.stats.multitest import multipletests

    valid = np.isfinite(raw)
    expected = np.full(raw.shape, np.nan)
    expected[valid] = multipletests(raw[valid], method="fdr_bh")[1]
    np.testing.assert_allclose(
        result["adjusted_p_value"],
        expected,
        equal_nan=True,
    )
    assert result["metadata"]["valid_count"] == int(valid.sum())
