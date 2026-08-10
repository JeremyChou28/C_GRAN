from __future__ import annotations

import json
import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, Optional

import numpy as np


@dataclass(frozen=True)
class PairChunk:
    """A contiguous chunk in the legacy lower-triangle pair order."""

    start: int
    stop: int
    row_i: np.ndarray
    row_j: np.ndarray
    packed_indices: np.ndarray


@dataclass(frozen=True)
class ArraySpec:
    path: Path
    dtype: np.dtype
    shape: tuple[int, ...]


def packed_upper_size(n_rows: int) -> int:
    if not isinstance(n_rows, int) or n_rows < 0:
        raise ValueError("n_rows must be a non-negative integer")
    return n_rows * (n_rows + 1) // 2


def pair_count(n_rows: int) -> int:
    if not isinstance(n_rows, int) or n_rows < 0:
        raise ValueError("n_rows must be a non-negative integer")
    return n_rows * (n_rows - 1) // 2


def packed_upper_index(i: int, j: int, n_rows: int) -> int:
    """Index of (i, j) in row-major packed upper triangle, diagonal included."""
    if not (0 <= i < n_rows and 0 <= j < n_rows):
        raise IndexError("matrix index out of range")
    if i > j:
        i, j = j, i
    row_offset = i * n_rows - i * (i - 1) // 2
    return row_offset + (j - i)


def packed_upper_indices(
    row_i: np.ndarray,
    row_j: np.ndarray,
    n_rows: int,
) -> np.ndarray:
    """Vectorized packed-upper indices for arbitrary symmetric positions."""
    row_i = np.asarray(row_i, dtype=np.int64)
    row_j = np.asarray(row_j, dtype=np.int64)
    if row_i.shape != row_j.shape:
        raise ValueError("row_i and row_j must have the same shape")
    lo = np.minimum(row_i, row_j)
    hi = np.maximum(row_i, row_j)
    if np.any(lo < 0) or np.any(hi >= n_rows):
        raise IndexError("matrix index out of range")
    return lo * n_rows - lo * (lo - 1) // 2 + (hi - lo)


def packed_symmetric_get(
    packed: np.ndarray,
    i: int,
    j: int,
    n_rows: int,
) -> Any:
    return packed[packed_upper_index(i, j, n_rows)]


def unpack_upper_triangle(
    packed: np.ndarray,
    n_rows: int,
    *,
    dtype: Optional[Any] = None,
) -> np.ndarray:
    expected = packed_upper_size(n_rows)
    if packed.ndim != 1 or packed.size != expected:
        raise ValueError("packed array has the wrong size")
    out = np.empty(
        (n_rows, n_rows),
        dtype=packed.dtype if dtype is None else dtype,
    )
    cursor = 0
    for i in range(n_rows):
        length = n_rows - i
        values = packed[cursor : cursor + length]
        out[i, i:] = values
        out[i:, i] = values
        cursor += length
    return out


def iter_legacy_pair_chunks(
    n_rows: int,
    chunk_size: int,
) -> Iterator[PairChunk]:
    """
    Yield non-diagonal pairs in the old script's order:

        for i in range(n_rows):
            for j in range(i):
                yield i, j

    Each chunk also contains indices into the engine's packed upper triangle.
    """
    if not isinstance(chunk_size, int) or chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    total = pair_count(n_rows)
    if total == 0:
        return

    row_i = np.empty(min(chunk_size, total), dtype=np.int64)
    row_j = np.empty(min(chunk_size, total), dtype=np.int64)
    filled = 0
    output_start = 0

    for i in range(1, n_rows):
        j_start = 0
        while j_start < i:
            capacity = row_i.size - filled
            take = min(capacity, i - j_start)
            sl = slice(filled, filled + take)
            row_i[sl] = i
            row_j[sl] = np.arange(j_start, j_start + take, dtype=np.int64)
            filled += take
            j_start += take

            if filled == row_i.size:
                i_chunk = row_i.copy()
                j_chunk = row_j.copy()
                yield PairChunk(
                    start=output_start,
                    stop=output_start + filled,
                    row_i=i_chunk,
                    row_j=j_chunk,
                    packed_indices=packed_upper_indices(
                        i_chunk,
                        j_chunk,
                        n_rows,
                    ),
                )
                output_start += filled
                remaining = total - output_start
                if remaining == 0:
                    return
                new_size = min(chunk_size, remaining)
                if new_size != row_i.size:
                    row_i = np.empty(new_size, dtype=np.int64)
                    row_j = np.empty(new_size, dtype=np.int64)
                filled = 0

    if filled:
        i_chunk = row_i[:filled].copy()
        j_chunk = row_j[:filled].copy()
        yield PairChunk(
            start=output_start,
            stop=output_start + filled,
            row_i=i_chunk,
            row_j=j_chunk,
            packed_indices=packed_upper_indices(
                i_chunk,
                j_chunk,
                n_rows,
            ),
        )


def create_memmap(
    path: os.PathLike[str] | str,
    *,
    dtype: Any,
    shape: tuple[int, ...],
    overwrite: bool = False,
) -> np.memmap:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not overwrite:
            raise FileExistsError(f"output already exists: {path}")
        path.unlink()
    return np.memmap(path, dtype=np.dtype(dtype), mode="w+", shape=shape)


def open_memmap(spec: ArraySpec, *, mode: str = "r") -> np.memmap:
    return np.memmap(spec.path, dtype=spec.dtype, mode=mode, shape=spec.shape)


def write_json_atomic(path: os.PathLike[str] | str, payload: Dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    temporary.replace(path)


def read_json(path: os.PathLike[str] | str) -> Dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as f:
        return json.load(f)


def check_free_space(
    directory: os.PathLike[str] | str,
    required_bytes: int,
    *,
    margin: float = 0.05,
) -> None:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(directory).free
    required_with_margin = math.ceil(required_bytes * (1.0 + margin))
    if free_bytes < required_with_margin:
        raise OSError(
            "insufficient disk space: "
            f"required_with_margin={required_with_margin}, "
            f"free={free_bytes}, directory={directory}"
        )


def flush_memmaps(*arrays: np.ndarray) -> None:
    for array in arrays:
        if isinstance(array, np.memmap):
            array.flush()


def close_memmaps(*arrays: np.ndarray) -> None:
    """Flush and close numpy memmap objects when possible."""
    for array in arrays:
        if isinstance(array, np.memmap):
            array.flush()
            mmap_obj = getattr(array, "_mmap", None)
            if mmap_obj is not None:
                mmap_obj.close()
