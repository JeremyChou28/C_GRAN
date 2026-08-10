from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from io_utils import (
    derive_before_bh_path,
    filter_outliers_iqr,
    load_intensity_table,
    sanitize_prefix,
    write_pair_csv,
)
from storage import pair_count
from tests.reference_legacy import filter_out_exceptions


def test_numpy_iqr_filter_matches_original_pandas_logic():
    rng = np.random.default_rng(21)
    values = rng.normal(size=(20, 17))
    values[rng.random(values.shape) < 0.2] = np.nan
    values[0] = np.nan
    values[1] = 4.0
    values[2, :] = np.nan
    values[2, :2] = [1.0, 100.0]
    values[3, -1] = 1000.0

    expected = filter_out_exceptions(pd.DataFrame(values.copy())).to_numpy()
    actual = filter_outliers_iqr(values)
    np.testing.assert_allclose(actual, expected, equal_nan=True)


def test_load_intensity_table_preserves_string_ids(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text(
        "ID\tS1\tS2\n 28773_2 \t1\t2\nABC\t3\t\n",
        encoding="utf-8",
    )
    loaded = load_intensity_table(
        path,
        expected_rows=2,
        expected_columns=2,
    )
    assert loaded.substance_ids.tolist() == ["28773_2", "ABC"]
    assert loaded.values.shape == (2, 2)
    assert np.isnan(loaded.values[1, 1])


def test_load_rejects_duplicate_or_empty_ids(tmp_path):
    duplicate = tmp_path / "duplicate.tsv"
    duplicate.write_text("ID\tS1\nA\t1\nA\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="重复"):
        load_intensity_table(duplicate)

    empty = tmp_path / "empty.tsv"
    empty.write_text("ID\tS1\n \t1\nB\t2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="空物质 ID"):
        load_intensity_table(empty)


def test_write_pair_csv_full_and_adjusted_filter(tmp_path):
    ids = np.array(["A", "B", "C", "D"])
    total = pair_count(len(ids))
    corr = np.linspace(-0.6, 0.6, total)
    p = np.linspace(0.01, 0.2, total)
    adjusted = np.array([0.01, 0.2, 0.03, np.nan, 0.049, 0.05])

    before = tmp_path / "before.csv"
    final = tmp_path / "final.csv"
    written_before = write_pair_csv(
        substance_ids=ids,
        correlations=corr,
        p_values=p,
        output_path=before,
        pair_chunk_size=2,
        show_progress=False,
    )
    written_final = write_pair_csv(
        substance_ids=ids,
        correlations=corr,
        p_values=p,
        adjusted_p_values=adjusted,
        adjusted_threshold=0.05,
        output_path=final,
        pair_chunk_size=2,
        show_progress=False,
    )

    before_df = pd.read_csv(before, dtype={"Substance 1": str, "Substance 2": str})
    final_df = pd.read_csv(final, dtype={"Substance 1": str, "Substance 2": str})
    assert written_before == total
    assert written_final == 3
    assert list(zip(before_df["Substance 1"], before_df["Substance 2"])) == [
        ("B", "A"),
        ("C", "A"),
        ("C", "B"),
        ("D", "A"),
        ("D", "B"),
        ("D", "C"),
    ]
    assert final_df["adjusted p-value"].tolist() == pytest.approx(
        [0.01, 0.03, 0.049]
    )


def test_path_and_prefix_helpers():
    assert derive_before_bh_path("out/result.csv").name == "result_before_bh.csv"
    assert sanitize_prefix("coral/MSV0001 pos") == "coral_MSV0001_pos"
