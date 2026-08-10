from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.stats import norm
from tqdm import tqdm

from storage import write_json_atomic


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

    Fisher/z calculations are performed only for pairs whose filtered and
    unfiltered correlations differ under the original nine-decimal/NaN rules.
    """
    filtered_corr = np.asarray(filtered_corr, dtype=np.float64)
    filtered_p = np.asarray(filtered_p, dtype=np.float64)
    filtered_count = np.asarray(filtered_count, dtype=np.float64)
    unfiltered_corr = np.asarray(unfiltered_corr, dtype=np.float64)
    unfiltered_p = np.asarray(unfiltered_p, dtype=np.float64)
    unfiltered_count = np.asarray(unfiltered_count, dtype=np.float64)

    if not (
        filtered_corr.shape
        == filtered_p.shape
        == filtered_count.shape
        == unfiltered_corr.shape
        == unfiltered_p.shape
        == unfiltered_count.shape
    ):
        raise ValueError("all merge inputs must have the same shape")

    final_corr = unfiltered_corr.copy()
    final_p = unfiltered_p.copy()

    both_missing = np.isnan(filtered_corr) & np.isnan(unfiltered_corr)
    one_missing = np.isnan(filtered_corr) | np.isnan(unfiltered_corr)
    values_different = (
        ~one_missing
        & (np.round(filtered_corr, 9) != np.round(unfiltered_corr, 9))
    )
    changed = (~both_missing) & (one_missing | values_different)
    changed_flat = np.flatnonzero(changed.ravel())
    if changed_flat.size == 0:
        return final_corr, final_p

    shape = filtered_corr.shape
    fc = filtered_corr.ravel()[changed_flat]
    fp = filtered_p.ravel()[changed_flat]
    fn = filtered_count.ravel()[changed_flat]
    uc = unfiltered_corr.ravel()[changed_flat]
    up = unfiltered_p.ravel()[changed_flat]
    un = unfiltered_count.ravel()[changed_flat]

    with np.errstate(divide="ignore", invalid="ignore"):
        z1 = np.where(
            fn <= 3,
            np.nan,
            legacy_fisher_z_transform(fc, fn),
        )
        z2 = np.where(
            un <= 3,
            np.nan,
            legacy_fisher_z_transform(uc, un),
        )
        difference_p = np.where(
            (fn <= 3) | (un <= 3),
            np.nan,
            legacy_z_test(z1, z2, fn, un),
        )

    replace = (difference_p < 0.05) | np.isnan(difference_p)
    if np.any(replace):
        target = changed_flat[replace]
        fc = fc[replace]
        fp = fp[replace]
        uc = uc[replace]
        up = up[replace]

        use_filtered = np.abs(fc) < np.abs(uc)
        selected_corr = np.where(use_filtered, fc, uc)
        selected_p = np.where(use_filtered, fp, up)

        # Preserve pandas Series.fillna(unfiltered) behavior from the old code.
        corr_flat = final_corr.ravel()
        p_flat = final_p.ravel()
        corr_valid = ~np.isnan(selected_corr)
        p_valid = ~np.isnan(selected_p)
        corr_flat[target[corr_valid]] = selected_corr[corr_valid]
        p_flat[target[p_valid]] = selected_p[p_valid]

    return final_corr.reshape(shape), final_p.reshape(shape)


def adjust_pvalues_bh(
    raw_p_values: np.ndarray,
    *,
    alpha: float = 0.05,
    scan_chunk_size: int = 10_000_000,
    metadata_path: Optional[os.PathLike[str] | str] = None,
    show_progress: bool = True,
) -> Dict[str, Any]:
    """
    Compute exact global BH results needed for ``adjusted p < alpha`` while
    using memory proportional to the number of BH-rejected hypotheses.

    The ordinary statsmodels implementation sorts every finite p-value and
    creates several arrays of that full length. For billions of tests this can
    require hundreds of GiB. This implementation instead:

    1. scans the raw p-values to count finite values;
    2. finds the exact BH rejection count by the monotone fixed-point rule
       ``R <- count(p <= alpha * R / m)``;
    3. extracts and sorts only the ``R`` rejected p-values;
    4. computes exact BH-adjusted values for that sorted significant prefix.

    The returned representation is sufficient to reproduce exactly every
    adjusted p-value that can be <= ``alpha``. It intentionally does not
    materialize adjusted values for hypotheses that cannot enter the final
    significant CSV.
    """
    raw_p_values = np.asarray(raw_p_values)
    if raw_p_values.ndim != 1:
        raise ValueError("raw_p_values must be one-dimensional")
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be between 0 and 1")
    if not isinstance(scan_chunk_size, int) or scan_chunk_size <= 0:
        raise ValueError("scan_chunk_size must be a positive integer")

    total_count = int(raw_p_values.size)
    metadata: Dict[str, Any] = {
        "complete": False,
        "method": "fdr_bh",
        "algorithm": "exact_significant_prefix_fixed_point",
        "alpha": float(alpha),
        "total_count": total_count,
        "scan_chunk_size": int(scan_chunk_size),
        "dtype": "float64",
    }
    metadata_path_str: Optional[str] = None
    if metadata_path is not None:
        path = Path(metadata_path)
        write_json_atomic(path, metadata)
        metadata_path_str = str(path.resolve())

    valid_count, below_alpha = _scan_valid_and_leq(
        raw_p_values,
        cutoff=alpha,
        chunk_size=scan_chunk_size,
        description="BH scan: valid and p<=alpha",
        show_progress=show_progress,
    )

    iteration_records = [
        {
            "iteration": 0,
            "input_count": int(valid_count),
            "cutoff": float(alpha),
            "output_count": int(below_alpha),
        }
    ]

    if valid_count == 0 or below_alpha == 0:
        rejection_count = 0
        cutoff = 0.0
    else:
        current = int(below_alpha)
        iteration = 1
        while True:
            cutoff = alpha * current / valid_count
            new_count = _scan_count_leq(
                raw_p_values,
                cutoff=cutoff,
                chunk_size=scan_chunk_size,
                description=(
                    f"BH fixed point {iteration}: "
                    f"p<={cutoff:.3e}"
                ),
                show_progress=show_progress,
            )
            iteration_records.append(
                {
                    "iteration": iteration,
                    "input_count": int(current),
                    "cutoff": float(cutoff),
                    "output_count": int(new_count),
                }
            )
            if new_count > current:
                raise RuntimeError(
                    "BH fixed-point count increased unexpectedly: "
                    f"current={current}, new={new_count}"
                )
            if new_count == current:
                break
            current = int(new_count)
            iteration += 1
            if current == 0:
                cutoff = 0.0
                break
        rejection_count = int(current)
        cutoff = (
            alpha * rejection_count / valid_count
            if rejection_count
            else 0.0
        )

    sorted_significant = _collect_values_leq(
        raw_p_values,
        cutoff=cutoff,
        expected_count=rejection_count,
        chunk_size=scan_chunk_size,
        show_progress=show_progress,
    )
    if sorted_significant.size:
        sorted_significant.sort(kind="quicksort")
        adjusted_significant = _adjust_sorted_bh_prefix(
            sorted_significant,
            valid_count=valid_count,
            chunk_size=scan_chunk_size,
        )
    else:
        adjusted_significant = np.empty(0, dtype=np.float64)

    metadata.update(
        {
            "complete": True,
            "valid_count": int(valid_count),
            "rejection_count": int(rejection_count),
            "raw_cutoff": float(cutoff),
            "fixed_point_iterations": iteration_records,
            "representation": (
                "sorted_significant_pvalues_plus_sorted_adjusted"
            ),
        }
    )
    if metadata_path is not None:
        write_json_atomic(Path(metadata_path), metadata)

    return {
        "sorted_significant_pvalues": sorted_significant,
        "adjusted_significant": adjusted_significant,
        "valid_count": int(valid_count),
        "rejection_count": int(rejection_count),
        "raw_cutoff": float(cutoff),
        "alpha": float(alpha),
        "metadata": metadata,
        "metadata_path": metadata_path_str,
    }


def lookup_bh_adjusted(
    p_values: np.ndarray,
    sorted_significant_pvalues: np.ndarray,
    adjusted_significant: np.ndarray,
) -> np.ndarray:
    """Look up exact BH-adjusted values for p-values in the rejected prefix."""
    p_values = np.asarray(p_values, dtype=np.float64)
    sorted_significant_pvalues = np.asarray(
        sorted_significant_pvalues,
        dtype=np.float64,
    )
    adjusted_significant = np.asarray(adjusted_significant, dtype=np.float64)
    if sorted_significant_pvalues.ndim != 1:
        raise ValueError("sorted_significant_pvalues must be one-dimensional")
    if adjusted_significant.ndim != 1:
        raise ValueError("adjusted_significant must be one-dimensional")
    if sorted_significant_pvalues.size != adjusted_significant.size:
        raise ValueError("BH lookup arrays must have the same size")
    if p_values.size == 0:
        return np.empty(p_values.shape, dtype=np.float64)
    if sorted_significant_pvalues.size == 0:
        raise ValueError("cannot look up non-empty p-values in an empty BH result")

    indices = np.searchsorted(
        sorted_significant_pvalues,
        p_values,
        side="right",
    ) - 1
    if np.any(indices < 0):
        raise ValueError("a p-value is below the stored significant range")

    # Values originate from the same raw p-value memmap, so exact equality is
    # expected. This check catches accidental lookup of a non-rejected value.
    matched = sorted_significant_pvalues[indices]
    if not np.array_equal(matched, p_values):
        raise ValueError("a p-value is not present in the BH significant prefix")
    return adjusted_significant[indices]


def _iter_slices(size: int, chunk_size: int):
    for start in range(0, size, chunk_size):
        yield start, min(start + chunk_size, size)


def _scan_valid_and_leq(
    values: np.ndarray,
    *,
    cutoff: float,
    chunk_size: int,
    description: str,
    show_progress: bool,
) -> Tuple[int, int]:
    valid_count = 0
    below_count = 0
    progress = tqdm(
        total=values.size,
        desc=description,
        disable=not show_progress,
    )
    try:
        for start, stop in _iter_slices(values.size, chunk_size):
            block = values[start:stop]
            valid_count += int(np.count_nonzero(np.isfinite(block)))
            # NaN compares False. Generated p-values never contain infinities.
            below_count += int(np.count_nonzero(block <= cutoff))
            progress.update(stop - start)
    finally:
        progress.close()
    return valid_count, below_count


def _scan_count_leq(
    values: np.ndarray,
    *,
    cutoff: float,
    chunk_size: int,
    description: str,
    show_progress: bool,
) -> int:
    count = 0
    progress = tqdm(
        total=values.size,
        desc=description,
        disable=not show_progress,
    )
    try:
        for start, stop in _iter_slices(values.size, chunk_size):
            count += int(np.count_nonzero(values[start:stop] <= cutoff))
            progress.update(stop - start)
    finally:
        progress.close()
    return count


def _collect_values_leq(
    values: np.ndarray,
    *,
    cutoff: float,
    expected_count: int,
    chunk_size: int,
    show_progress: bool,
) -> np.ndarray:
    output = np.empty(expected_count, dtype=np.float64)
    if expected_count == 0:
        return output

    cursor = 0
    progress = tqdm(
        total=values.size,
        desc="BH collect significant p-values",
        disable=not show_progress,
    )
    try:
        for start, stop in _iter_slices(values.size, chunk_size):
            block = values[start:stop]
            selected = np.asarray(block[block <= cutoff], dtype=np.float64)
            next_cursor = cursor + selected.size
            if next_cursor > expected_count:
                raise RuntimeError(
                    "collected more BH-significant p-values than expected"
                )
            output[cursor:next_cursor] = selected
            cursor = next_cursor
            progress.update(stop - start)
    finally:
        progress.close()

    if cursor != expected_count:
        raise RuntimeError(
            "BH significant count changed between scans: "
            f"expected={expected_count}, collected={cursor}"
        )
    return output


def _adjust_sorted_bh_prefix(
    sorted_pvalues: np.ndarray,
    *,
    valid_count: int,
    chunk_size: int,
) -> np.ndarray:
    count = sorted_pvalues.size
    adjusted = np.empty(count, dtype=np.float64)

    for start, stop in _iter_slices(count, chunk_size):
        ranks = np.arange(start + 1, stop + 1, dtype=np.float64)
        adjusted[start:stop] = (
            sorted_pvalues[start:stop] * float(valid_count) / ranks
        )

    running_min = np.inf
    stop = count
    while stop > 0:
        start = max(0, stop - chunk_size)
        block = adjusted[start:stop]
        reverse = block[::-1]
        np.minimum.accumulate(reverse, out=reverse)
        if np.isfinite(running_min):
            np.minimum(block, running_min, out=block)
        running_min = float(block[0])
        stop = start

    np.clip(adjusted, 0.0, 1.0, out=adjusted)
    return adjusted
