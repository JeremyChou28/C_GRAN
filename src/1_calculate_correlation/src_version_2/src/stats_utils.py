from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from storage import create_memmap, flush_memmaps, write_json_atomic


def legacy_fisher_z_transform(
    correlation: np.ndarray,
    common_count: np.ndarray,
) -> np.ndarray:
    """Vectorized form of the existing script's Fisher transformation."""
    correlation = np.asarray(correlation, dtype=np.float64)
    common_count = np.asarray(common_count, dtype=np.float64)
    valid_r = np.clip(correlation, -0.999999999, 0.999999999)
    return 0.5 * (
        0.5 * (np.log1p(valid_r) - np.log1p(-valid_r))
        / np.sqrt(1.0 / (common_count - 3.0))
    )


def legacy_z_test(
    z1: np.ndarray,
    z2: np.ndarray,
    n1: np.ndarray,
    n2: np.ndarray,
) -> np.ndarray:
    """Vectorized form of the existing script's z-test."""
    z1 = np.asarray(z1, dtype=np.float64)
    z2 = np.asarray(z2, dtype=np.float64)
    n1 = np.asarray(n1, dtype=np.float64)
    n2 = np.asarray(n2, dtype=np.float64)
    se = np.sqrt(1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0))
    z_statistic = (z1 - z2) / se
    return 2.0 * (1.0 - norm.cdf(np.abs(z_statistic)))


def merge_filtered_unfiltered_chunk(
    filtered_corr: np.ndarray,
    filtered_p: np.ndarray,
    filtered_count: np.ndarray,
    unfiltered_corr: np.ndarray,
    unfiltered_p: np.ndarray,
    unfiltered_count: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reproduce the old compare -> Fisher -> min-|r| replacement logic.

    The behavior for NaN, n <= 3, nine-decimal comparison and replacement is
    intentionally kept compatible with the original script.
    """
    filtered_corr = np.asarray(filtered_corr, dtype=np.float64)
    filtered_p = np.asarray(filtered_p, dtype=np.float64)
    filtered_count = np.asarray(filtered_count, dtype=np.float64)
    unfiltered_corr = np.asarray(unfiltered_corr, dtype=np.float64)
    unfiltered_p = np.asarray(unfiltered_p, dtype=np.float64)
    unfiltered_count = np.asarray(unfiltered_count, dtype=np.float64)

    both_missing = np.isnan(filtered_corr) & np.isnan(unfiltered_corr)
    one_missing = np.isnan(filtered_corr) | np.isnan(unfiltered_corr)
    values_different = (
        ~one_missing
        & (
            np.round(filtered_corr, 9)
            != np.round(unfiltered_corr, 9)
        )
    )
    changed = (~both_missing) & (one_missing | values_different)

    with np.errstate(divide="ignore", invalid="ignore"):
        z1 = np.where(
            filtered_count <= 3,
            np.nan,
            legacy_fisher_z_transform(filtered_corr, filtered_count),
        )
        z2 = np.where(
            unfiltered_count <= 3,
            np.nan,
            legacy_fisher_z_transform(unfiltered_corr, unfiltered_count),
        )
        difference_p = np.where(
            (filtered_count <= 3) | (unfiltered_count <= 3),
            np.nan,
            legacy_z_test(
                z1,
                z2,
                filtered_count,
                unfiltered_count,
            ),
        )

    significant_or_nan = changed & (
        (difference_p < 0.05) | np.isnan(difference_p)
    )

    use_filtered = np.abs(filtered_corr) < np.abs(unfiltered_corr)
    correlation_min = np.where(
        use_filtered,
        filtered_corr,
        unfiltered_corr,
    )
    p_value_min = np.where(
        use_filtered,
        filtered_p,
        unfiltered_p,
    )

    candidate_corr = np.where(significant_or_nan, correlation_min, np.nan)
    candidate_p = np.where(significant_or_nan, p_value_min, np.nan)

    # pandas Series.fillna(base) from the old implementation.
    final_corr = np.where(
        np.isnan(candidate_corr),
        unfiltered_corr,
        candidate_corr,
    )
    final_p = np.where(
        np.isnan(candidate_p),
        unfiltered_p,
        candidate_p,
    )
    return final_corr, final_p


def adjust_pvalues_bh(
    raw_p_values: np.ndarray,
    output_path: os.PathLike[str] | str,
    *,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """Run one global in-memory Benjamini-Hochberg correction."""
    raw_p_values = np.asarray(raw_p_values)
    if raw_p_values.ndim != 1:
        raise ValueError("raw_p_values must be one-dimensional")

    output_path = Path(output_path)
    adjusted = create_memmap(
        output_path,
        dtype=np.float64,
        shape=raw_p_values.shape,
        overwrite=overwrite,
    )
    adjusted[:] = np.nan

    valid_mask = np.isfinite(raw_p_values)
    valid_pvalues = np.asarray(raw_p_values[valid_mask], dtype=np.float64)
    if valid_pvalues.size:
        adjusted_valid = multipletests(
            valid_pvalues,
            method="fdr_bh",
        )[1]
        adjusted[valid_mask] = adjusted_valid
        del adjusted_valid
    flush_memmaps(adjusted)

    metadata = {
        "complete": True,
        "method": "fdr_bh",
        "total_count": int(raw_p_values.size),
        "valid_count": int(valid_pvalues.size),
        "path": str(output_path.resolve()),
        "dtype": "float64",
    }
    metadata_path = output_path.with_suffix(output_path.suffix + ".json")
    write_json_atomic(metadata_path, metadata)
    del valid_pvalues

    return {
        "adjusted_p_value": adjusted,
        "metadata": metadata,
        "metadata_path": str(metadata_path),
    }
