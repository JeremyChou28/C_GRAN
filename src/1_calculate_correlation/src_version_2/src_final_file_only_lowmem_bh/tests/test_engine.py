import warnings

import numpy as np
import pytest
from scipy.stats import pearsonr

from engine import calculate_correlations, open_result
from storage import unpack_upper_triangle


def reference(values: np.ndarray):
    n = values.shape[0]
    corr = np.full((n, n), np.nan, dtype=np.float64)
    p_value = np.full((n, n), np.nan, dtype=np.float64)
    count = np.zeros((n, n), dtype=np.uint64)
    for i in range(n):
        for j in range(i, n):
            common = ~np.isnan(values[i]) & ~np.isnan(values[j])
            c = int(common.sum())
            count[i, j] = count[j, i] = c
            if c < 3:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = pearsonr(values[i, common], values[j, common])
            corr[i, j] = corr[j, i] = result.statistic
            p_value[i, j] = p_value[j, i] = result.pvalue
    return corr, p_value, count


def assert_engine_matches(values, actual, atol=5e-10):
    expected_corr, expected_p, expected_count = reference(values)
    n = values.shape[0]
    if actual["metadata"]["triangular"]:
        corr = unpack_upper_triangle(actual["corr"], n)
        p_value = unpack_upper_triangle(actual["p_value"], n)
        count = unpack_upper_triangle(actual["common_count"], n)
    else:
        corr = np.asarray(actual["corr"])
        p_value = np.asarray(actual["p_value"])
        count = np.asarray(actual["common_count"])

    np.testing.assert_array_equal(count, expected_count)
    np.testing.assert_allclose(
        corr,
        expected_corr,
        atol=atol,
        rtol=atol,
        equal_nan=True,
    )

    endpoint = (
        np.isfinite(corr)
        & np.isfinite(expected_corr)
        & (np.abs(corr) >= 1 - 1e-12)
        & (np.abs(expected_corr) >= 1 - 1e-12)
    )
    np.testing.assert_allclose(
        p_value[~endpoint],
        expected_p[~endpoint],
        atol=1e-10,
        rtol=1e-7,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        p_value[endpoint],
        expected_p[endpoint],
        atol=2e-7,
        rtol=1e-7,
        equal_nan=True,
    )


def test_random_pairwise_complete_against_scipy():
    rng = np.random.default_rng(11)
    values = rng.normal(size=(36, 23))
    values[rng.random(values.shape) < 0.35] = np.nan
    values[0] = np.nan
    values[1, :6] = 3.0

    actual = calculate_correlations(
        values,
        jobs=3,
        storage="memory",
        block_size=7,
        triangular=True,
    )
    assert set(actual) == {
        "corr",
        "p_value",
        "common_count",
        "metadata",
        "metadata_path",
    }
    assert_engine_matches(values, actual)


def test_full_and_packed_outputs_are_equivalent():
    rng = np.random.default_rng(12)
    values = rng.normal(size=(20, 13))
    values[rng.random(values.shape) < 0.2] = np.nan

    packed = calculate_correlations(
        values,
        jobs=2,
        storage="memory",
        block_size=6,
        triangular=True,
    )
    full = calculate_correlations(
        values,
        jobs=2,
        storage="memory",
        block_size=6,
        triangular=False,
    )

    np.testing.assert_allclose(
        unpack_upper_triangle(packed["corr"], values.shape[0]),
        full["corr"],
        equal_nan=True,
    )
    np.testing.assert_allclose(
        unpack_upper_triangle(packed["p_value"], values.shape[0]),
        full["p_value"],
        equal_nan=True,
    )
    np.testing.assert_array_equal(
        unpack_upper_triangle(packed["common_count"], values.shape[0]),
        full["common_count"],
    )


def test_memmap_reopen_and_common_count_dtype(tmp_path):
    values = np.array(
        [
            [1.0, 2.0, np.nan, 4.0],
            [2.0, 4.0, 6.0, np.nan],
            [3.0, np.nan, 1.0, 8.0],
        ]
    )
    result = calculate_correlations(
        values,
        jobs=1,
        output_dir=tmp_path,
        prefix="case",
        block_size=2,
        overwrite=False,
    )
    assert isinstance(result["corr"], np.memmap)
    assert result["common_count"].dtype == np.uint8
    reopened = open_result(result["metadata_path"])
    assert_engine_matches(values, reopened)


def test_count_dtype_selection_and_validation():
    values = np.ones((2, 256), dtype=np.float64)
    result = calculate_correlations(
        values,
        jobs=1,
        storage="memory",
        triangular=True,
    )
    assert result["common_count"].dtype == np.uint16

    with pytest.raises(ValueError, match="cannot hold"):
        calculate_correlations(
            values,
            jobs=1,
            storage="memory",
            count_dtype=np.uint8,
        )


def test_constant_and_insufficient_pairs():
    values = np.array(
        [
            [1.0, 1.0, 1.0, np.nan],
            [2.0, 3.0, 4.0, 5.0],
            [2.0, np.nan, np.nan, 8.0],
        ]
    )
    result = calculate_correlations(
        values,
        jobs=1,
        storage="memory",
        triangular=False,
    )
    assert result["common_count"][0, 1] == 3
    assert np.isnan(result["corr"][0, 1])
    assert np.isnan(result["p_value"][0, 1])
    assert result["common_count"][0, 2] == 1
    assert np.isnan(result["corr"][0, 2])


def test_large_offset_small_variation_is_stable():
    base = 1e12
    values = np.array(
        [
            [base + 0.001, base + 0.002, base + 0.003, base + 0.004],
            [base + 0.002, base + 0.004, base + 0.006, base + 0.008],
        ]
    )
    result = calculate_correlations(
        values,
        jobs=1,
        storage="memory",
        triangular=False,
    )
    assert result["corr"][0, 1] > 0.999


def test_engine_rejects_dataframe_and_large_memory_output():
    import pandas as pd

    with pytest.raises(TypeError, match="numpy.ndarray"):
        calculate_correlations(
            pd.DataFrame([[1.0, 2.0, 3.0]]),
            jobs=1,
            storage="memory",
        )

    with pytest.raises(ValueError, match="N < 10_000"):
        calculate_correlations(
            np.zeros((10_000, 1), dtype=np.float64),
            jobs=1,
            storage="memory",
        )


def test_existing_memmap_files_are_not_overwritten(tmp_path):
    values = np.arange(12, dtype=np.float64).reshape(3, 4)
    calculate_correlations(
        values,
        jobs=1,
        output_dir=tmp_path,
        prefix="same",
    )
    with pytest.raises(FileExistsError):
        calculate_correlations(
            values,
            jobs=1,
            output_dir=tmp_path,
            prefix="same",
        )


def test_fused_engine_writes_legacy_pair_order_and_matches_reference(tmp_path):
    from engine import calculate_merged_correlations
    from io_utils import filter_outliers_iqr
    from tests.reference_legacy import run_legacy_pipeline

    rng = np.random.default_rng(1234)
    values = rng.normal(size=(9, 11))
    values[rng.random(values.shape) < 0.2] = np.nan
    values[3, -1] = 1000.0
    ids = np.array([f"ID_{i}" for i in range(values.shape[0])])
    expected_before, _, _ = run_legacy_pipeline(values, ids)

    result = calculate_merged_correlations(
        filter_outliers_iqr(values),
        values,
        jobs=3,
        block_size=4,
        output_dir=tmp_path,
        prefix="fused",
    )
    np.testing.assert_allclose(
        result["corr"],
        expected_before["Correlation"].to_numpy(),
        atol=5e-9,
        rtol=5e-9,
        equal_nan=True,
    )
    expected_p = expected_before["P-Value"].to_numpy()
    actual_p = np.asarray(result["p_value"])
    endpoint = (
        np.isfinite(result["corr"])
        & np.isfinite(expected_before["Correlation"].to_numpy())
        & (np.abs(result["corr"]) >= 1 - 1e-10)
        & (np.abs(expected_before["Correlation"].to_numpy()) >= 1 - 1e-10)
    )
    np.testing.assert_allclose(
        actual_p[~endpoint],
        expected_p[~endpoint],
        atol=1e-9,
        rtol=1e-7,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        actual_p[endpoint],
        expected_p[endpoint],
        atol=2e-7,
        rtol=1e-7,
        equal_nan=True,
    )
    assert result["metadata"]["pair_order"] == (
        "for i in range(n_rows): for j in range(i)"
    )
