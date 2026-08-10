from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from storage import iter_legacy_pair_chunks, pair_count


@dataclass(frozen=True)
class InputMatrix:
    substance_ids: np.ndarray
    values: np.ndarray


def normalize_substance_ids(values) -> np.ndarray:
    series = pd.Series(values, dtype="string").str.strip()
    series = series.mask(series == "")
    if series.isna().any():
        invalid_count = int(series.isna().sum())
        raise ValueError(
            f"输入文件第一列中存在 {invalid_count} 个空物质 ID。"
        )
    if series.duplicated().any():
        examples = series[series.duplicated()].drop_duplicates().head(10).tolist()
        raise ValueError(f"输入文件中存在重复物质 ID。示例：{examples}")
    return series.astype(str).to_numpy()


def load_intensity_table(
    path: os.PathLike[str] | str,
    *,
    expected_rows: Optional[int] = None,
    expected_columns: Optional[int] = None,
) -> InputMatrix:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    substance_ids = normalize_substance_ids(frame.index)

    if expected_rows is not None and frame.shape[0] != expected_rows:
        raise ValueError(
            "The shape read from the file should be "
            "(number of substances, number of samples). "
            f"实际 shape={frame.shape}，compounds_num={expected_rows}, "
            f"samples_num={expected_columns}。"
        )
    if expected_columns is not None and frame.shape[1] != expected_columns:
        raise ValueError(
            "The shape read from the file should be "
            "(number of substances, number of samples). "
            f"实际 shape={frame.shape}，compounds_num={expected_rows}, "
            f"samples_num={expected_columns}。"
        )

    try:
        numeric = frame.apply(pd.to_numeric, errors="raise")
    except Exception as exc:
        raise ValueError("丰度区域包含无法转换为数值的内容。") from exc

    values = numeric.to_numpy(dtype=np.float64, copy=True)
    if np.isinf(values).any():
        raise ValueError("丰度矩阵包含正无穷或负无穷；缺失值必须使用 NaN。")
    return InputMatrix(substance_ids=substance_ids, values=values)



def load_substance_ids(
    path: os.PathLike[str] | str,
    *,
    expected_rows: Optional[int] = None,
) -> np.ndarray:
    """Load and validate only the first (substance-ID) column of a TSV."""
    frame = pd.read_csv(path, sep="\t", usecols=[0], dtype="string")
    substance_ids = normalize_substance_ids(frame.iloc[:, 0])
    if expected_rows is not None and substance_ids.size != expected_rows:
        raise ValueError(
            f"expected {expected_rows} substance IDs, got {substance_ids.size}"
        )
    return substance_ids

def filter_outliers_iqr(mat: np.ndarray) -> np.ndarray:
    """NumPy equivalent of the original row-wise pandas IQR filtering."""
    if not isinstance(mat, np.ndarray) or mat.ndim != 2:
        raise ValueError("mat must be a two-dimensional numpy array")
    values = np.asarray(mat, dtype=np.float64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        q1 = np.nanquantile(values, 0.25, axis=1, method="linear")
        q3 = np.nanquantile(values, 0.75, axis=1, method="linear")
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    filtered = values.copy()
    outlier_mask = (
        (filtered < lower[:, None])
        | (filtered > upper[:, None])
    )
    filtered[outlier_mask] = np.nan
    return filtered


def sanitize_prefix(value: str) -> str:
    value = str(value).strip()
    if not value:
        raise ValueError("tmp_name must not be empty")
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return sanitized.strip("._-") or "result"


def derive_before_bh_path(final_path: os.PathLike[str] | str) -> Path:
    final_path = Path(final_path)
    suffix = final_path.suffix or ".csv"
    stem = final_path.stem if final_path.suffix else final_path.name
    return final_path.with_name(f"{stem}_before_bh{suffix}")


def write_pair_csv(
    *,
    substance_ids: np.ndarray,
    correlations: np.ndarray,
    p_values: np.ndarray,
    output_path: os.PathLike[str] | str,
    pair_chunk_size: int,
    adjusted_p_values: Optional[np.ndarray] = None,
    adjusted_threshold: Optional[float] = None,
    overwrite: bool = False,
    show_progress: bool = True,
) -> int:
    """
    Write pair arrays in the legacy lower-triangle order using DataFrame chunks.

    When ``adjusted_threshold`` is provided, only rows with a finite adjusted
    p-value below that threshold are written.
    """
    substance_ids = np.asarray(substance_ids)
    n_rows = substance_ids.size
    expected_pairs = pair_count(n_rows)
    correlations = np.asarray(correlations)
    p_values = np.asarray(p_values)
    if correlations.ndim != 1 or correlations.size != expected_pairs:
        raise ValueError("correlations has the wrong shape")
    if p_values.ndim != 1 or p_values.size != expected_pairs:
        raise ValueError("p_values has the wrong shape")
    if adjusted_p_values is not None:
        adjusted_p_values = np.asarray(adjusted_p_values)
        if adjusted_p_values.ndim != 1 or adjusted_p_values.size != expected_pairs:
            raise ValueError("adjusted_p_values has the wrong shape")
    if adjusted_threshold is not None and adjusted_p_values is None:
        raise ValueError(
            "adjusted_p_values is required when adjusted_threshold is set"
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {output_path}")
    temporary = output_path.with_suffix(output_path.suffix + ".partial")
    if temporary.exists():
        if not overwrite:
            raise FileExistsError(f"partial output already exists: {temporary}")
        temporary.unlink()

    columns = ["Substance 1", "Substance 2", "Correlation", "P-Value"]
    if adjusted_p_values is not None:
        columns.append("adjusted p-value")

    written = 0
    first_chunk = True
    progress = tqdm(
        total=expected_pairs,
        desc=f"Writing {output_path.name}",
        disable=not show_progress,
    )
    try:
        for chunk in iter_legacy_pair_chunks(n_rows, pair_chunk_size):
            corr = correlations[chunk.start : chunk.stop]
            p_value = p_values[chunk.start : chunk.stop]
            adjusted = (
                None
                if adjusted_p_values is None
                else adjusted_p_values[chunk.start : chunk.stop]
            )

            if adjusted_threshold is None:
                selection = slice(None)
                selected_count = chunk.stop - chunk.start
            else:
                selection = np.isfinite(adjusted) & (
                    adjusted < adjusted_threshold
                )
                selected_count = int(np.count_nonzero(selection))

            if selected_count:
                data = {
                    "Substance 1": substance_ids[chunk.row_i[selection]],
                    "Substance 2": substance_ids[chunk.row_j[selection]],
                    "Correlation": corr[selection],
                    "P-Value": p_value[selection],
                }
                if adjusted is not None:
                    data["adjusted p-value"] = adjusted[selection]
                pd.DataFrame(data, columns=columns).to_csv(
                    temporary,
                    mode="w" if first_chunk else "a",
                    header=first_chunk,
                    index=False,
                )
                first_chunk = False
                written += selected_count
            progress.update(chunk.stop - chunk.start)

        if first_chunk:
            pd.DataFrame(columns=columns).to_csv(temporary, index=False)
        temporary.replace(output_path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    finally:
        progress.close()
    return written


def write_bh_significant_csv(
    *,
    substance_ids: np.ndarray,
    correlations: np.ndarray,
    p_values: np.ndarray,
    sorted_significant_pvalues: np.ndarray,
    adjusted_significant: np.ndarray,
    raw_cutoff: float,
    rejection_count: int,
    adjusted_threshold: float,
    output_path: os.PathLike[str] | str,
    pair_chunk_size: int,
    overwrite: bool = False,
    show_progress: bool = True,
) -> int:
    """
    Write pairs whose exact global BH-adjusted p-value is below a threshold.

    Only the BH-rejected p-value prefix is stored in memory. During this
    sequential pass, raw p-values at or below ``raw_cutoff`` are mapped to
    their exact adjusted values through the sorted significant prefix.
    """
    from stats_utils import lookup_bh_adjusted

    substance_ids = np.asarray(substance_ids)
    n_rows = substance_ids.size
    expected_pairs = pair_count(n_rows)
    correlations = np.asarray(correlations)
    p_values = np.asarray(p_values)
    sorted_significant_pvalues = np.asarray(
        sorted_significant_pvalues,
        dtype=np.float64,
    )
    adjusted_significant = np.asarray(
        adjusted_significant,
        dtype=np.float64,
    )

    if correlations.ndim != 1 or correlations.size != expected_pairs:
        raise ValueError("correlations has the wrong shape")
    if p_values.ndim != 1 or p_values.size != expected_pairs:
        raise ValueError("p_values has the wrong shape")
    if sorted_significant_pvalues.ndim != 1:
        raise ValueError("sorted_significant_pvalues must be one-dimensional")
    if adjusted_significant.ndim != 1:
        raise ValueError("adjusted_significant must be one-dimensional")
    if sorted_significant_pvalues.size != adjusted_significant.size:
        raise ValueError("BH arrays must have the same size")
    if sorted_significant_pvalues.size != rejection_count:
        raise ValueError("rejection_count does not match the BH arrays")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {output_path}")
    temporary = output_path.with_suffix(output_path.suffix + ".partial")
    if temporary.exists():
        if not overwrite:
            raise FileExistsError(f"partial output already exists: {temporary}")
        temporary.unlink()

    columns = [
        "Substance 1",
        "Substance 2",
        "Correlation",
        "P-Value",
        "adjusted p-value",
    ]
    candidate_count = 0
    written = 0
    first_chunk = True
    progress = tqdm(
        total=expected_pairs,
        desc=f"Writing {output_path.name}",
        disable=not show_progress,
    )
    try:
        for chunk in iter_legacy_pair_chunks(n_rows, pair_chunk_size):
            chunk_p = p_values[chunk.start : chunk.stop]
            candidate = chunk_p <= raw_cutoff
            selected_positions = np.flatnonzero(candidate)
            candidate_count += int(selected_positions.size)

            if selected_positions.size:
                selected_p = np.asarray(
                    chunk_p[selected_positions],
                    dtype=np.float64,
                )
                selected_adjusted = lookup_bh_adjusted(
                    selected_p,
                    sorted_significant_pvalues,
                    adjusted_significant,
                )
                output_selection = np.isfinite(selected_adjusted) & (
                    selected_adjusted < adjusted_threshold
                )
                output_positions = selected_positions[output_selection]
                selected_count = int(output_positions.size)
            else:
                selected_adjusted = np.empty(0, dtype=np.float64)
                output_selection = np.empty(0, dtype=bool)
                output_positions = np.empty(0, dtype=np.int64)
                selected_count = 0

            if selected_count:
                data = {
                    "Substance 1": substance_ids[
                        chunk.row_i[output_positions]
                    ],
                    "Substance 2": substance_ids[
                        chunk.row_j[output_positions]
                    ],
                    "Correlation": correlations[
                        chunk.start : chunk.stop
                    ][output_positions],
                    "P-Value": chunk_p[output_positions],
                    "adjusted p-value": selected_adjusted[
                        output_selection
                    ],
                }
                pd.DataFrame(data, columns=columns).to_csv(
                    temporary,
                    mode="w" if first_chunk else "a",
                    header=first_chunk,
                    index=False,
                )
                first_chunk = False
                written += selected_count
            progress.update(chunk.stop - chunk.start)

        if candidate_count != rejection_count:
            raise RuntimeError(
                "BH candidate count changed between scans: "
                f"expected={rejection_count}, observed={candidate_count}"
            )
        if first_chunk:
            pd.DataFrame(columns=columns).to_csv(temporary, index=False)
        temporary.replace(output_path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    finally:
        progress.close()
    return written
