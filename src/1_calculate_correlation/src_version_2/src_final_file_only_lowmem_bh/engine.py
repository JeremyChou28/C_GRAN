from __future__ import annotations

import math
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterator, Optional, Tuple

import numpy as np
from scipy.special import betainc
from threadpoolctl import threadpool_limits

from stats_utils import merge_filtered_unfiltered_chunk

from storage import (
    check_free_space,
    create_memmap,
    flush_memmaps,
    packed_upper_index,
    packed_upper_size,
    read_json,
    write_json_atomic,
)


_RESULT_VERSION = 2
_MEMORY_ROW_LIMIT = 10_000


def calculate_correlations(
    mat: np.ndarray,
    jobs: int,
    backend: str = "cpu",
    *,
    block_size: Optional[int] = None,
    output_dtype: Any = np.float64,
    count_dtype: Optional[Any] = None,
    storage: str = "memmap",
    output_dir: Optional[os.PathLike[str] | str] = None,
    prefix: str = "result",
    triangular: bool = True,
    blas_threads: Optional[int] = None,
    check_disk_space: bool = True,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """
    Calculate pairwise-complete Pearson correlations between matrix rows.

    Only features that are non-NaN in both rows are used. Intermediate
    calculations are float64. The result contains ``corr``, ``p_value`` and
    ``common_count``. Correlation and p-value are NaN when fewer than three
    common observations exist or when either row is constant on the pairwise
    intersection.

    By default, arrays are disk-backed memmaps containing a packed upper
    triangle, diagonal included. In-memory output is only supported for
    ``N < 10_000``.
    """
    mat = _validate_inputs(
        mat=mat,
        jobs=jobs,
        backend=backend,
        output_dtype=output_dtype,
        count_dtype=count_dtype,
        storage=storage,
        block_size=block_size,
        prefix=prefix,
    )
    n_rows, n_features = mat.shape
    output_dtype = np.dtype(output_dtype)
    count_dtype = _resolve_count_dtype(n_features, count_dtype)
    storage = storage.lower()

    if block_size is None:
        block_size = _choose_block_size(n_rows=n_rows, jobs=jobs)
    block_size = min(int(block_size), n_rows)

    tasks = list(_iter_upper_block_tasks(n_rows, block_size))
    worker_count = min(jobs, max(1, len(tasks)))

    cpu_count = os.cpu_count() or 1
    if blas_threads is None:
        blas_threads = max(1, cpu_count // worker_count)
    if not isinstance(blas_threads, int) or blas_threads <= 0:
        raise ValueError("blas_threads must be a positive integer")

    values, values_sq, mask_f = _prepare_data(mat)

    outputs, metadata, metadata_path = _allocate_outputs(
        n_rows=n_rows,
        n_features=n_features,
        output_dtype=output_dtype,
        count_dtype=count_dtype,
        storage=storage,
        output_dir=output_dir,
        prefix=prefix,
        triangular=triangular,
        block_size=block_size,
        jobs=jobs,
        blas_threads=blas_threads,
        check_disk_space_enabled=check_disk_space,
        overwrite=overwrite,
    )

    corr_out = outputs["corr"]
    p_out = outputs["p_value"]
    count_out = outputs["common_count"]

    try:
        with threadpool_limits(limits=blas_threads, user_api="blas"):
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = [
                    executor.submit(
                        _compute_and_store_block,
                        task,
                        values,
                        values_sq,
                        mask_f,
                        corr_out,
                        p_out,
                        count_out,
                        n_rows,
                        triangular,
                    )
                    for task in tasks
                ]
                for future in as_completed(futures):
                    future.result()

        flush_memmaps(corr_out, p_out, count_out)
        metadata["complete"] = True
        if metadata_path is not None:
            write_json_atomic(metadata_path, metadata)
    except Exception:
        metadata["complete"] = False
        if metadata_path is not None:
            write_json_atomic(metadata_path, metadata)
        raise

    return {
        "corr": corr_out,
        "p_value": p_out,
        "common_count": count_out,
        "metadata": metadata,
        "metadata_path": metadata_path,
    }


def estimate_output_bytes(
    n_rows: int,
    output_dtype: Any = np.float64,
    *,
    count_dtype: Any = np.uint16,
    triangular: bool = True,
) -> int:
    if not isinstance(n_rows, int) or n_rows < 0:
        raise ValueError("n_rows must be a non-negative integer")
    output_dtype = np.dtype(output_dtype)
    count_dtype = np.dtype(count_dtype)
    if output_dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
        raise ValueError("output_dtype must be float32 or float64")
    if not np.issubdtype(count_dtype, np.unsignedinteger):
        raise ValueError("count_dtype must be an unsigned integer dtype")
    elements = packed_upper_size(n_rows) if triangular else n_rows * n_rows
    return elements * (2 * output_dtype.itemsize + count_dtype.itemsize)


def open_result(
    metadata_path: os.PathLike[str] | str,
    *,
    mode: str = "r",
) -> Dict[str, Any]:
    metadata = read_json(metadata_path)
    if metadata.get("storage") != "memmap":
        raise ValueError("metadata does not describe a memmap result")
    if not metadata.get("complete", False):
        raise ValueError("result is incomplete")

    n_rows = int(metadata["n_rows"])
    triangular = bool(metadata["triangular"])
    shape = (
        (packed_upper_size(n_rows),)
        if triangular
        else (n_rows, n_rows)
    )
    paths = metadata["paths"]
    return {
        "corr": np.memmap(
            paths["corr"],
            dtype=np.dtype(metadata["output_dtype"]),
            mode=mode,
            shape=shape,
        ),
        "p_value": np.memmap(
            paths["p_value"],
            dtype=np.dtype(metadata["output_dtype"]),
            mode=mode,
            shape=shape,
        ),
        "common_count": np.memmap(
            paths["common_count"],
            dtype=np.dtype(metadata["count_dtype"]),
            mode=mode,
            shape=shape,
        ),
        "metadata": metadata,
        "metadata_path": str(metadata_path),
    }


def _validate_inputs(
    *,
    mat: np.ndarray,
    jobs: int,
    backend: str,
    output_dtype: Any,
    count_dtype: Optional[Any],
    storage: str,
    block_size: Optional[int],
    prefix: str,
) -> np.ndarray:
    if not isinstance(mat, np.ndarray):
        raise TypeError("mat must be a numpy.ndarray")
    if mat.ndim != 2:
        raise ValueError("mat must be a two-dimensional array")
    if mat.shape[0] <= 0:
        raise ValueError("mat must contain at least one row")
    if mat.shape[1] <= 0:
        raise ValueError("mat must contain at least one feature")
    if not np.issubdtype(mat.dtype, np.number):
        raise TypeError("mat must contain numeric values")
    if np.isinf(mat).any():
        raise ValueError(
            "mat contains positive or negative infinity; only np.nan is "
            "accepted as a missing value"
        )

    if not isinstance(jobs, int) or isinstance(jobs, bool) or jobs <= 0:
        raise ValueError("jobs must be a positive integer")
    if not isinstance(backend, str) or backend.lower() != "cpu":
        raise ValueError("backend must be 'cpu'")
    if not isinstance(storage, str) or storage.lower() not in {
        "memory",
        "memmap",
    }:
        raise ValueError("storage must be 'memory' or 'memmap'")
    if storage.lower() == "memory" and mat.shape[0] >= _MEMORY_ROW_LIMIT:
        raise ValueError(
            "storage='memory' is only supported when N < 10_000; "
            "use storage='memmap'"
        )

    dtype = np.dtype(output_dtype)
    if dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
        raise ValueError("output_dtype must be float32 or float64")
    if count_dtype is not None:
        _resolve_count_dtype(mat.shape[1], count_dtype)
    if block_size is not None and (
        not isinstance(block_size, int)
        or isinstance(block_size, bool)
        or block_size <= 0
    ):
        raise ValueError("block_size must be a positive integer")
    if not isinstance(prefix, str) or not prefix.strip():
        raise ValueError("prefix must be a non-empty string")
    return mat


def _resolve_count_dtype(n_features: int, count_dtype: Optional[Any]) -> np.dtype:
    if count_dtype is None:
        if n_features <= np.iinfo(np.uint8).max:
            return np.dtype(np.uint8)
        if n_features <= np.iinfo(np.uint16).max:
            return np.dtype(np.uint16)
        if n_features <= np.iinfo(np.uint32).max:
            return np.dtype(np.uint32)
        return np.dtype(np.uint64)

    dtype = np.dtype(count_dtype)
    if not np.issubdtype(dtype, np.unsignedinteger):
        raise ValueError("count_dtype must be an unsigned integer dtype")
    if n_features > np.iinfo(dtype).max:
        raise ValueError(
            f"count_dtype={dtype.name} cannot hold n_features={n_features}"
        )
    return dtype


def _choose_block_size(n_rows: int, jobs: int) -> int:
    total_workspace_budget = 2 * 1024**3
    per_worker_budget = min(
        192 * 1024**2,
        max(32 * 1024**2, total_workspace_budget // max(1, jobs)),
    )
    raw = int(math.sqrt(per_worker_budget / 96.0))
    candidates = (256, 512, 1024, 2048)
    chosen = max(c for c in candidates if c <= max(256, raw))
    return min(chosen, n_rows)


def _prepare_data(
    mat: np.ndarray,
    *,
    row_chunk_size: int = 4096,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_rows, n_features = mat.shape
    values = np.empty((n_rows, n_features), dtype=np.float64)
    mask_f = np.empty((n_rows, n_features), dtype=np.float64)

    for start in range(0, n_rows, row_chunk_size):
        stop = min(start + row_chunk_size, n_rows)
        block = np.asarray(mat[start:stop], dtype=np.float64)
        valid = ~np.isnan(block)
        counts = valid.sum(axis=1)

        first_index = np.argmax(valid, axis=1)
        row_index = np.arange(stop - start)
        pivots = block[row_index, first_index].copy()
        pivots[counts == 0] = 0.0

        centered = block - pivots[:, None]
        centered[~valid] = 0.0

        mean_delta = np.zeros(stop - start, dtype=np.float64)
        np.divide(
            centered.sum(axis=1, dtype=np.float64),
            counts,
            out=mean_delta,
            where=counts > 0,
        )
        centered -= mean_delta[:, None]
        centered[~valid] = 0.0

        scales = np.max(np.abs(centered), axis=1)
        np.divide(
            centered,
            scales[:, None],
            out=centered,
            where=scales[:, None] > 0,
        )
        centered[~valid] = 0.0

        if not np.isfinite(centered).all():
            raise FloatingPointError(
                "row normalization produced a non-finite value"
            )

        values[start:stop] = centered
        mask_f[start:stop] = valid

    return values, np.square(values), mask_f


def _iter_upper_block_tasks(
    n_rows: int,
    block_size: int,
) -> Iterator[Tuple[int, int, int, int]]:
    for i0 in range(0, n_rows, block_size):
        i1 = min(i0 + block_size, n_rows)
        for j0 in range(i0, n_rows, block_size):
            j1 = min(j0 + block_size, n_rows)
            yield i0, i1, j0, j1


def _compute_block(
    values_i: np.ndarray,
    values_j: np.ndarray,
    values_sq_i: np.ndarray,
    values_sq_j: np.ndarray,
    mask_i: np.ndarray,
    mask_j: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    count = mask_i @ mask_j.T
    sum_x = values_i @ mask_j.T
    sum_y = mask_i @ values_j.T
    sum_x2 = values_sq_i @ mask_j.T
    sum_y2 = mask_i @ values_sq_j.T
    sum_xy = values_i @ values_j.T

    illegal = count < 3.0
    inv_count = np.zeros_like(count)
    np.divide(1.0, count, out=inv_count, where=count > 0)

    temp = np.empty_like(count)
    np.multiply(sum_x, sum_y, out=temp)
    temp *= inv_count
    sum_xy -= temp

    np.multiply(sum_x, sum_x, out=temp)
    temp *= inv_count
    sum_x2 -= temp

    np.multiply(sum_y, sum_y, out=temp)
    temp *= inv_count
    sum_y2 -= temp

    zero_tol = 64.0 * np.finfo(np.float64).eps * np.maximum(count, 1.0)
    sum_x2[sum_x2 <= zero_tol] = 0.0
    sum_y2[sum_y2 <= zero_tol] = 0.0

    np.multiply(sum_x2, sum_y2, out=temp)
    np.sqrt(temp, out=temp)
    valid_corr = (~illegal) & (temp > 0.0)

    corr = sum_xy
    corr[~valid_corr] = np.nan
    np.divide(corr, temp, out=corr, where=valid_corr)
    np.clip(corr, -1.0, 1.0, out=corr)

    p_value = np.full_like(corr, np.nan)
    if np.any(valid_corr):
        r_valid = corr[valid_corr]
        beta_x = np.maximum(1.0 - r_valid * r_valid, 0.0)
        beta_a = (count[valid_corr] - 2.0) / 2.0
        p_value[valid_corr] = betainc(beta_a, 0.5, beta_x)

    return corr, p_value, count


def _compute_and_store_block(
    task: Tuple[int, int, int, int],
    values: np.ndarray,
    values_sq: np.ndarray,
    mask_f: np.ndarray,
    corr_out: np.ndarray,
    p_out: np.ndarray,
    count_out: np.ndarray,
    n_rows: int,
    triangular: bool,
) -> None:
    i0, i1, j0, j1 = task
    corr, p_value, count = _compute_block(
        values[i0:i1],
        values[j0:j1],
        values_sq[i0:i1],
        values_sq[j0:j1],
        mask_f[i0:i1],
        mask_f[j0:j1],
    )

    if i0 == j0:
        _force_symmetric_tile(corr)
        _force_symmetric_tile(p_value)
        _force_symmetric_tile(count)

    if triangular:
        _write_packed_upper_block(corr_out, corr, i0, i1, j0, j1, n_rows)
        _write_packed_upper_block(p_out, p_value, i0, i1, j0, j1, n_rows)
        _write_packed_upper_block(count_out, count, i0, i1, j0, j1, n_rows)
    else:
        corr_out[i0:i1, j0:j1] = corr
        p_out[i0:i1, j0:j1] = p_value
        count_out[i0:i1, j0:j1] = count
        if i0 != j0:
            corr_out[j0:j1, i0:i1] = corr.T
            p_out[j0:j1, i0:i1] = p_value.T
            count_out[j0:j1, i0:i1] = count.T


def _force_symmetric_tile(tile: np.ndarray) -> None:
    for row in range(1, tile.shape[0]):
        tile[row, :row] = tile[:row, row]


def _write_packed_upper_block(
    destination: np.ndarray,
    block: np.ndarray,
    i0: int,
    i1: int,
    j0: int,
    j1: int,
    n_rows: int,
) -> None:
    if i0 < j0:
        for local_i, global_i in enumerate(range(i0, i1)):
            start = packed_upper_index(global_i, j0, n_rows)
            destination[start : start + (j1 - j0)] = block[local_i]
    else:
        for local_i, global_i in enumerate(range(i0, i1)):
            local_start = local_i
            global_j = j0 + local_start
            start = packed_upper_index(global_i, global_j, n_rows)
            values = block[local_i, local_start:]
            destination[start : start + values.size] = values


def _allocate_outputs(
    *,
    n_rows: int,
    n_features: int,
    output_dtype: np.dtype,
    count_dtype: np.dtype,
    storage: str,
    output_dir: Optional[os.PathLike[str] | str],
    prefix: str,
    triangular: bool,
    block_size: int,
    jobs: int,
    blas_threads: int,
    check_disk_space_enabled: bool,
    overwrite: bool,
) -> Tuple[Dict[str, np.ndarray], Dict[str, Any], Optional[str]]:
    shape = (
        (packed_upper_size(n_rows),)
        if triangular
        else (n_rows, n_rows)
    )
    metadata: Dict[str, Any] = {
        "version": _RESULT_VERSION,
        "complete": False,
        "storage": storage,
        "triangular": triangular,
        "n_rows": n_rows,
        "n_features": n_features,
        "matrix_shape": [n_rows, n_rows],
        "storage_shape": list(shape),
        "output_dtype": output_dtype.name,
        "count_dtype": count_dtype.name,
        "block_size": block_size,
        "jobs": jobs,
        "blas_threads": blas_threads,
        "prefix": prefix,
        "paths": None,
    }

    if storage == "memory":
        return (
            {
                "corr": np.empty(shape, dtype=output_dtype),
                "p_value": np.empty(shape, dtype=output_dtype),
                "common_count": np.empty(shape, dtype=count_dtype),
            },
            metadata,
            None,
        )

    directory = Path.cwd() if output_dir is None else Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {
        "corr": str((directory / f"{prefix}_corr.dat").resolve()),
        "p_value": str((directory / f"{prefix}_p_value.dat").resolve()),
        "common_count": str(
            (directory / f"{prefix}_common_count.dat").resolve()
        ),
    }
    metadata_path = str((directory / f"{prefix}_metadata.json").resolve())
    metadata["paths"] = paths
    metadata["metadata_path"] = metadata_path
    required_bytes = estimate_output_bytes(
        n_rows,
        output_dtype,
        count_dtype=count_dtype,
        triangular=triangular,
    )
    metadata["estimated_output_bytes"] = required_bytes

    existing = [Path(path) for path in paths.values()] + [Path(metadata_path)]
    if not overwrite:
        collision = next((path for path in existing if path.exists()), None)
        if collision is not None:
            raise FileExistsError(f"output already exists: {collision}")
    else:
        for path in existing:
            if path.exists():
                path.unlink()

    if check_disk_space_enabled:
        check_free_space(directory, required_bytes)

    outputs = {
        "corr": create_memmap(
            paths["corr"], dtype=output_dtype, shape=shape, overwrite=False
        ),
        "p_value": create_memmap(
            paths["p_value"], dtype=output_dtype, shape=shape, overwrite=False
        ),
        "common_count": create_memmap(
            paths["common_count"], dtype=count_dtype, shape=shape, overwrite=False
        ),
    }
    write_json_atomic(metadata_path, metadata)
    return outputs, metadata, metadata_path



def calculate_merged_correlations(
    filtered_mat: np.ndarray,
    unfiltered_mat: np.ndarray,
    jobs: int,
    backend: str = "cpu",
    *,
    block_size: Optional[int] = None,
    output_dtype: Any = np.float64,
    output_dir: os.PathLike[str] | str,
    prefix: str = "merged",
    blas_threads: Optional[int] = None,
    check_disk_space: bool = True,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """
    Calculate filtered and unfiltered correlations block-by-block and merge
    them immediately with the legacy comparison/Fisher/min-|r| logic.

    Unlike calling :func:`calculate_correlations` twice, this function never
    materializes filtered or unfiltered N x N result arrays. Only two
    non-diagonal pair-order memmaps are written:

    ``corr`` and ``p_value`` in the order::

        for i in range(n_rows):
            for j in range(i):
                yield (i, j)

    This removes the separate O(N^2) "Merging pairs" scan and its six large
    input memmaps. ``output_dtype`` controls the quantization applied to each
    filtered/unfiltered block before the legacy merge, matching the previous
    two-engine pipeline. The merged outputs remain float64, as before.
    """
    filtered_mat = _validate_inputs(
        mat=filtered_mat,
        jobs=jobs,
        backend=backend,
        output_dtype=output_dtype,
        count_dtype=None,
        storage="memmap",
        block_size=block_size,
        prefix=prefix,
    )
    unfiltered_mat = _validate_inputs(
        mat=unfiltered_mat,
        jobs=jobs,
        backend=backend,
        output_dtype=output_dtype,
        count_dtype=None,
        storage="memmap",
        block_size=block_size,
        prefix=prefix,
    )
    if filtered_mat.shape != unfiltered_mat.shape:
        raise ValueError(
            "filtered_mat and unfiltered_mat must have the same shape"
        )

    n_rows, n_features = unfiltered_mat.shape
    storage_dtype = np.dtype(output_dtype)
    if block_size is None:
        block_size = _choose_block_size(n_rows=n_rows, jobs=jobs)
    block_size = min(int(block_size), n_rows)

    tasks = list(_iter_upper_block_tasks(n_rows, block_size))
    worker_count = min(jobs, max(1, len(tasks)))
    cpu_count = os.cpu_count() or 1
    if blas_threads is None:
        blas_threads = max(1, cpu_count // worker_count)
    if not isinstance(blas_threads, int) or blas_threads <= 0:
        raise ValueError("blas_threads must be a positive integer")

    # Row normalization is cheap relative to O(N^2 M), and keeps the same
    # stable sufficient-statistics computation used by calculate_correlations.
    filtered_values, filtered_sq, filtered_mask = _prepare_data(filtered_mat)
    unfiltered_values, unfiltered_sq, unfiltered_mask = _prepare_data(
        unfiltered_mat
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pair_total = n_rows * (n_rows - 1) // 2
    paths = {
        "corr": output_dir / f"{prefix}_corr.dat",
        "p_value": output_dir / f"{prefix}_p_value.dat",
        "metadata": output_dir / f"{prefix}_metadata.json",
    }
    existing = list(paths.values())
    if overwrite:
        for path in existing:
            if path.exists():
                path.unlink()
    else:
        collision = next((path for path in existing if path.exists()), None)
        if collision is not None:
            raise FileExistsError(f"output already exists: {collision}")

    required_bytes = pair_total * 2 * np.dtype(np.float64).itemsize
    if check_disk_space:
        check_free_space(output_dir, required_bytes)

    corr_out = create_memmap(
        paths["corr"], dtype=np.float64, shape=(pair_total,), overwrite=False
    )
    p_out = create_memmap(
        paths["p_value"],
        dtype=np.float64,
        shape=(pair_total,),
        overwrite=False,
    )
    metadata: Dict[str, Any] = {
        "version": _RESULT_VERSION,
        "complete": False,
        "kind": "fused_filtered_unfiltered_legacy_merge",
        "n_rows": n_rows,
        "n_features": n_features,
        "pair_count": pair_total,
        "pair_order": "for i in range(n_rows): for j in range(i)",
        "source_output_dtype": storage_dtype.name,
        "merged_output_dtype": "float64",
        "block_size": block_size,
        "jobs": jobs,
        "blas_threads": blas_threads,
        "prefix": prefix,
        "estimated_output_bytes": required_bytes,
        "paths": {key: str(path.resolve()) for key, path in paths.items()},
    }
    write_json_atomic(paths["metadata"], metadata)

    try:
        with threadpool_limits(limits=blas_threads, user_api="blas"):
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = [
                    executor.submit(
                        _compute_merge_and_store_block,
                        task,
                        filtered_values,
                        filtered_sq,
                        filtered_mask,
                        unfiltered_values,
                        unfiltered_sq,
                        unfiltered_mask,
                        corr_out,
                        p_out,
                        storage_dtype,
                    )
                    for task in tasks
                ]
                for future in as_completed(futures):
                    future.result()

        flush_memmaps(corr_out, p_out)
        metadata["complete"] = True
        write_json_atomic(paths["metadata"], metadata)
    except Exception:
        metadata["complete"] = False
        write_json_atomic(paths["metadata"], metadata)
        raise

    return {
        "corr": corr_out,
        "p_value": p_out,
        "metadata": metadata,
        "metadata_path": str(paths["metadata"].resolve()),
    }


def _compute_merge_and_store_block(
    task: Tuple[int, int, int, int],
    filtered_values: np.ndarray,
    filtered_sq: np.ndarray,
    filtered_mask: np.ndarray,
    unfiltered_values: np.ndarray,
    unfiltered_sq: np.ndarray,
    unfiltered_mask: np.ndarray,
    corr_out: np.ndarray,
    p_out: np.ndarray,
    source_output_dtype: np.dtype,
) -> None:
    i0, i1, j0, j1 = task
    filtered_corr, filtered_p, filtered_count = _compute_block(
        filtered_values[i0:i1],
        filtered_values[j0:j1],
        filtered_sq[i0:i1],
        filtered_sq[j0:j1],
        filtered_mask[i0:i1],
        filtered_mask[j0:j1],
    )
    unfiltered_corr, unfiltered_p, unfiltered_count = _compute_block(
        unfiltered_values[i0:i1],
        unfiltered_values[j0:j1],
        unfiltered_sq[i0:i1],
        unfiltered_sq[j0:j1],
        unfiltered_mask[i0:i1],
        unfiltered_mask[j0:j1],
    )

    # The previous pipeline stored each source result in output_dtype before
    # reading it back for merge. Reproduce that quantization exactly.
    if source_output_dtype == np.dtype(np.float32):
        filtered_corr = filtered_corr.astype(np.float32).astype(np.float64)
        filtered_p = filtered_p.astype(np.float32).astype(np.float64)
        unfiltered_corr = unfiltered_corr.astype(np.float32).astype(np.float64)
        unfiltered_p = unfiltered_p.astype(np.float32).astype(np.float64)

    merged_corr, merged_p = merge_filtered_unfiltered_chunk(
        filtered_corr,
        filtered_p,
        filtered_count,
        unfiltered_corr,
        unfiltered_p,
        unfiltered_count,
    )
    _write_legacy_pairs_from_upper_block(
        corr_out,
        p_out,
        merged_corr,
        merged_p,
        i0,
        i1,
        j0,
        j1,
    )


def _write_legacy_pairs_from_upper_block(
    corr_out: np.ndarray,
    p_out: np.ndarray,
    corr_block: np.ndarray,
    p_block: np.ndarray,
    i0: int,
    i1: int,
    j0: int,
    j1: int,
) -> None:
    """Write one upper block into contiguous legacy lower-pair segments."""
    if i0 < j0:
        width = i1 - i0
        for local_j, global_j in enumerate(range(j0, j1)):
            start = global_j * (global_j - 1) // 2 + i0
            corr_out[start : start + width] = corr_block[:, local_j]
            p_out[start : start + width] = p_block[:, local_j]
        return

    # Diagonal block: for each higher row j, store lower columns i0 <= i < j.
    for local_j in range(1, j1 - j0):
        global_j = j0 + local_j
        start = global_j * (global_j - 1) // 2 + i0
        corr_out[start : start + local_j] = corr_block[:local_j, local_j]
        p_out[start : start + local_j] = p_block[:local_j, local_j]
