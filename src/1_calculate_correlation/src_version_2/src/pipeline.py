from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from tqdm import tqdm

from engine import calculate_correlations
from io_utils import (
    derive_before_bh_path,
    filter_outliers_iqr,
    load_intensity_table,
    sanitize_prefix,
    write_pair_csv,
)
from stats_utils import adjust_pvalues_bh, merge_filtered_unfiltered_chunk
from storage import (
    create_memmap,
    flush_memmaps,
    iter_legacy_pair_chunks,
    pair_count,
    write_json_atomic,
)


@dataclass(frozen=True)
class PipelineConfig:
    intensity_file: Path
    compounds_num: int
    samples_num: int
    correlation_result_filename: Path
    n_jobs: int
    tmp_name: str
    block_size: Optional[int] = None
    output_dtype: Any = np.float64
    blas_threads: Optional[int] = None
    pair_chunk_size: int = 1_000_000
    save_before_bh_csv: bool = True
    overwrite: bool = False
    show_progress: bool = True


@dataclass(frozen=True)
class PipelineResult:
    before_bh_csv: Optional[Path]
    significant_csv: Path
    filtered_metadata: Path
    unfiltered_metadata: Path
    merged_metadata: Path
    bh_metadata: Path
    total_pairs: int
    significant_pairs: int


def run_pipeline(config: PipelineConfig) -> PipelineResult:
    _validate_config(config)

    final_csv = Path(config.correlation_result_filename).resolve()
    before_bh_csv = (
        derive_before_bh_path(final_csv)
        if config.save_before_bh_csv
        else None
    )
    output_dir = final_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = sanitize_prefix(config.tmp_name)

    _check_csv_collisions(
        before_bh_csv=before_bh_csv,
        final_csv=final_csv,
        overwrite=config.overwrite,
    )

    print(f"Loading data {config.intensity_file}...")
    input_matrix = load_intensity_table(
        config.intensity_file,
        expected_rows=config.compounds_num,
        expected_columns=config.samples_num,
    )
    substance_ids = input_matrix.substance_ids
    values = input_matrix.values
    n_rows = values.shape[0]
    total_pairs = pair_count(n_rows)
    print(f"Data loaded. shape: {values.shape}")
    print("Substance ID examples:", substance_ids[:5].tolist())

    print("Filtering row-wise IQR outliers...")
    filtered_values = filter_outliers_iqr(values)

    print("Calculating filtered correlations...")
    filtered_result = calculate_correlations(
        filtered_values,
        jobs=config.n_jobs,
        backend="cpu",
        block_size=config.block_size,
        output_dtype=config.output_dtype,
        storage="memmap",
        output_dir=output_dir,
        prefix=f"{prefix}_filtered",
        triangular=True,
        blas_threads=config.blas_threads,
        overwrite=config.overwrite,
    )
    del filtered_values

    print("Calculating unfiltered correlations...")
    unfiltered_result = calculate_correlations(
        values,
        jobs=config.n_jobs,
        backend="cpu",
        block_size=config.block_size,
        output_dtype=config.output_dtype,
        storage="memmap",
        output_dir=output_dir,
        prefix=f"{prefix}_unfiltered",
        triangular=True,
        blas_threads=config.blas_threads,
        overwrite=config.overwrite,
    )
    del values

    merged_paths = {
        "corr": output_dir / f"{prefix}_merged_corr.dat",
        "p_value": output_dir / f"{prefix}_merged_p_value.dat",
        "adjusted_p_value": output_dir
        / f"{prefix}_merged_adjusted_p_value.dat",
        "metadata": output_dir / f"{prefix}_merged_metadata.json",
    }
    _prepare_merged_paths(merged_paths, overwrite=config.overwrite)

    merged_corr = create_memmap(
        merged_paths["corr"],
        dtype=np.float64,
        shape=(total_pairs,),
        overwrite=False,
    )
    merged_p = create_memmap(
        merged_paths["p_value"],
        dtype=np.float64,
        shape=(total_pairs,),
        overwrite=False,
    )

    merged_metadata: Dict[str, Any] = {
        "complete": False,
        "n_rows": n_rows,
        "pair_count": total_pairs,
        "pair_order": "for i in range(n_rows): for j in range(i)",
        "dtype": "float64",
        "paths": {
            key: str(path.resolve())
            for key, path in merged_paths.items()
        },
    }
    write_json_atomic(merged_paths["metadata"], merged_metadata)

    print("Merging filtered and unfiltered results...")
    progress = tqdm(
        total=total_pairs,
        desc="Merging pairs",
        disable=not config.show_progress,
    )
    try:
        for chunk in iter_legacy_pair_chunks(
            n_rows,
            config.pair_chunk_size,
        ):
            packed = chunk.packed_indices
            final_corr, final_p = merge_filtered_unfiltered_chunk(
                filtered_result["corr"][packed],
                filtered_result["p_value"][packed],
                filtered_result["common_count"][packed],
                unfiltered_result["corr"][packed],
                unfiltered_result["p_value"][packed],
                unfiltered_result["common_count"][packed],
            )
            merged_corr[chunk.start : chunk.stop] = final_corr
            merged_p[chunk.start : chunk.stop] = final_p
            progress.update(chunk.stop - chunk.start)
        flush_memmaps(merged_corr, merged_p)
        merged_metadata["complete"] = True
        write_json_atomic(merged_paths["metadata"], merged_metadata)
    except Exception:
        merged_metadata["complete"] = False
        write_json_atomic(merged_paths["metadata"], merged_metadata)
        raise
    finally:
        progress.close()

    if before_bh_csv is not None:
        print(f"Writing pre-BH merged CSV: {before_bh_csv}")
        write_pair_csv(
            substance_ids=substance_ids,
            correlations=merged_corr,
            p_values=merged_p,
            output_path=before_bh_csv,
            pair_chunk_size=config.pair_chunk_size,
            overwrite=config.overwrite,
            show_progress=config.show_progress,
        )
    else:
        print("Skipping pre-BH merged CSV.")

    print("Running global Benjamini-Hochberg correction...")
    bh_result = adjust_pvalues_bh(
        merged_p,
        merged_paths["adjusted_p_value"],
        overwrite=False,
    )

    print(f"Writing adjusted-p < 0.05 CSV: {final_csv}")
    significant_pairs = write_pair_csv(
        substance_ids=substance_ids,
        correlations=merged_corr,
        p_values=merged_p,
        adjusted_p_values=bh_result["adjusted_p_value"],
        adjusted_threshold=0.05,
        output_path=final_csv,
        pair_chunk_size=config.pair_chunk_size,
        overwrite=config.overwrite,
        show_progress=config.show_progress,
    )

    return PipelineResult(
        before_bh_csv=before_bh_csv,
        significant_csv=final_csv,
        filtered_metadata=Path(filtered_result["metadata_path"]),
        unfiltered_metadata=Path(unfiltered_result["metadata_path"]),
        merged_metadata=merged_paths["metadata"],
        bh_metadata=Path(bh_result["metadata_path"]),
        total_pairs=total_pairs,
        significant_pairs=significant_pairs,
    )


def _validate_config(config: PipelineConfig) -> None:
    if config.compounds_num <= 0:
        raise ValueError("compounds_num must be positive")
    if config.samples_num <= 0:
        raise ValueError("samples_num must be positive")
    if config.n_jobs <= 0:
        raise ValueError("n_jobs must be positive")
    if config.pair_chunk_size <= 0:
        raise ValueError("pair_chunk_size must be positive")
    dtype = np.dtype(config.output_dtype)
    if dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
        raise ValueError("output_dtype must be float32 or float64")


def _check_csv_collisions(
    before_bh_csv: Optional[Path],
    final_csv: Path,
    *,
    overwrite: bool,
) -> None:
    paths = [final_csv]
    if before_bh_csv is not None:
        paths.insert(0, before_bh_csv)

    for path in paths:
        partial = path.with_suffix(path.suffix + ".partial")
        if not overwrite:
            if path.exists():
                raise FileExistsError(f"output already exists: {path}")
            if partial.exists():
                raise FileExistsError(f"partial output already exists: {partial}")
        else:
            if path.exists():
                path.unlink()
            if partial.exists():
                partial.unlink()


def _prepare_merged_paths(
    paths: Dict[str, Path],
    *,
    overwrite: bool,
) -> None:
    for path in paths.values():
        sidecar = (
            path.with_suffix(path.suffix + ".json")
            if path.name.endswith("adjusted_p_value.dat")
            else None
        )
        candidates = [path] + ([] if sidecar is None else [sidecar])
        for candidate in candidates:
            if candidate.exists():
                if not overwrite:
                    raise FileExistsError(
                        f"output already exists: {candidate}"
                    )
                candidate.unlink()
