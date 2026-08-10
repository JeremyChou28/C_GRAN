from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import PipelineConfig, run_pipeline
from tests.reference_legacy import run_legacy_pipeline


def make_values():
    base = np.arange(1.0, 13.0)
    rng = np.random.default_rng(33)
    values = np.vstack(
        [
            base,
            2.0 * base,
            -base,
            rng.normal(size=12),
            base.copy(),
            rng.normal(size=12),
            base + rng.normal(scale=0.05, size=12),
            np.full(12, 4.0),
        ]
    )
    values[4, -1] = 1000.0
    values[5, [1, 5, 9]] = np.nan
    values[6, [0, 3]] = np.nan
    return values


def assert_csv_numeric_close(actual, expected, include_adjusted=False):
    np.testing.assert_array_equal(
        actual[["Substance 1", "Substance 2"]].astype(str).to_numpy(),
        expected[["Substance 1", "Substance 2"]].astype(str).to_numpy(),
    )
    np.testing.assert_allclose(
        actual["Correlation"],
        expected["Correlation"],
        atol=5e-9,
        rtol=5e-9,
        equal_nan=True,
    )
    endpoint = (
        np.isfinite(actual["Correlation"].to_numpy())
        & np.isfinite(expected["Correlation"].to_numpy())
        & (np.abs(actual["Correlation"].to_numpy()) >= 1 - 1e-10)
        & (np.abs(expected["Correlation"].to_numpy()) >= 1 - 1e-10)
    )
    np.testing.assert_allclose(
        actual.loc[~endpoint, "P-Value"],
        expected.loc[~endpoint, "P-Value"],
        atol=1e-9,
        rtol=1e-7,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        actual.loc[endpoint, "P-Value"],
        expected.loc[endpoint, "P-Value"],
        atol=2e-7,
        rtol=1e-7,
        equal_nan=True,
    )
    if include_adjusted:
        np.testing.assert_allclose(
            actual["adjusted p-value"],
            expected["adjusted p-value"],
            atol=2e-7,
            rtol=1e-6,
            equal_nan=True,
        )


def test_full_pipeline_matches_original_algorithm(tmp_path):
    values = make_values()
    ids = np.array([f"CMP_{i}" for i in range(values.shape[0])])
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=ids).to_csv(input_path, sep="\t")

    expected_before, expected_final, _ = run_legacy_pipeline(values, ids)

    final_path = tmp_path / "result.csv"
    result = run_pipeline(
        PipelineConfig(
            intensity_file=input_path,
            compounds_num=values.shape[0],
            samples_num=values.shape[1],
            correlation_result_filename=final_path,
            n_jobs=2,
            tmp_name="folder/test case",
            block_size=3,
            pair_chunk_size=5,
            overwrite=False,
            show_progress=False,
        )
    )

    actual_before = pd.read_csv(
        result.before_bh_csv,
        dtype={"Substance 1": str, "Substance 2": str},
    )
    actual_final = pd.read_csv(
        result.significant_csv,
        dtype={"Substance 1": str, "Substance 2": str},
    )

    assert result.before_bh_csv.name == "result_before_bh.csv"
    assert list(actual_before.columns) == [
        "Substance 1",
        "Substance 2",
        "Correlation",
        "P-Value",
    ]
    assert list(actual_final.columns) == [
        "Substance 1",
        "Substance 2",
        "Correlation",
        "P-Value",
        "adjusted p-value",
    ]
    assert_csv_numeric_close(actual_before, expected_before)
    assert_csv_numeric_close(actual_final, expected_final, include_adjusted=True)

    assert result.total_pairs == values.shape[0] * (values.shape[0] - 1) // 2
    assert result.significant_pairs == len(actual_final)
    assert result.filtered_metadata.exists()
    assert result.unfiltered_metadata.exists()
    assert result.merged_metadata.exists()
    assert result.bh_metadata.exists()

    output_names = {path.name for path in tmp_path.iterdir()}
    assert not any("illegal_mask" in name for name in output_names)
    assert any("common_count" in name for name in output_names)
    assert any("folder_test_case_filtered" in name for name in output_names)


def test_pipeline_refuses_to_overwrite_existing_outputs(tmp_path):
    values = make_values()[:4]
    ids = np.array([f"X_{i}" for i in range(values.shape[0])])
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=ids).to_csv(input_path, sep="\t")
    final_path = tmp_path / "result.csv"
    config = PipelineConfig(
        intensity_file=input_path,
        compounds_num=values.shape[0],
        samples_num=values.shape[1],
        correlation_result_filename=final_path,
        n_jobs=1,
        tmp_name="same",
        block_size=2,
        pair_chunk_size=3,
        show_progress=False,
    )
    run_pipeline(config)

    import pytest

    with pytest.raises(FileExistsError):
        run_pipeline(config)



def test_pipeline_can_skip_before_bh_csv(tmp_path):
    values = make_values()
    ids = np.array([f"SKIP_{i}" for i in range(values.shape[0])])
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=ids).to_csv(input_path, sep="\t")

    final_path = tmp_path / "result.csv"
    expected_before_path = tmp_path / "result_before_bh.csv"
    result = run_pipeline(
        PipelineConfig(
            intensity_file=input_path,
            compounds_num=values.shape[0],
            samples_num=values.shape[1],
            correlation_result_filename=final_path,
            n_jobs=2,
            tmp_name="skip-before-bh",
            block_size=3,
            pair_chunk_size=5,
            save_before_bh_csv=False,
            overwrite=False,
            show_progress=False,
        )
    )

    assert result.before_bh_csv is None
    assert not expected_before_path.exists()
    assert result.significant_csv == final_path.resolve()
    assert final_path.exists()
    assert result.bh_metadata.exists()


def test_skipped_before_bh_csv_does_not_participate_in_collision_check(tmp_path):
    values = make_values()[:4]
    ids = np.array([f"COLLISION_{i}" for i in range(values.shape[0])])
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=ids).to_csv(input_path, sep="\t")

    # A stale pre-BH file should not block a run that explicitly skips it.
    stale_before_bh = tmp_path / "result_before_bh.csv"
    stale_before_bh.write_text("stale\n", encoding="utf-8")

    final_path = tmp_path / "result.csv"
    result = run_pipeline(
        PipelineConfig(
            intensity_file=input_path,
            compounds_num=values.shape[0],
            samples_num=values.shape[1],
            correlation_result_filename=final_path,
            n_jobs=1,
            tmp_name="skip-collision",
            block_size=2,
            pair_chunk_size=3,
            save_before_bh_csv=False,
            show_progress=False,
        )
    )

    assert result.before_bh_csv is None
    assert stale_before_bh.read_text(encoding="utf-8") == "stale\n"
    assert final_path.exists()
