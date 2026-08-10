from __future__ import annotations

import gc
import shutil
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

from engine import calculate_merged_correlations
from io_utils import (
    derive_before_bh_path,
    filter_outliers_iqr,
    load_intensity_table,
    sanitize_prefix,
    write_bh_significant_csv,
    write_pair_csv,
)
from stats_utils import adjust_pvalues_bh
from storage import close_memmaps, pair_count


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
    bh_chunk_size: int = 10_000_000
    save_before_bh_csv: bool = True
    final_only: bool = False
    overwrite: bool = False
    show_progress: bool = True


@dataclass(frozen=True)
class PipelineResult:
    before_bh_csv: Optional[Path]
    significant_csv: Path
    work_dir: Optional[Path]
    merged_metadata: Optional[Path]
    bh_metadata: Optional[Path]
    total_pairs: int
    significant_pairs: int


def run_pipeline(config: PipelineConfig) -> PipelineResult:
    _validate_config(config)

    final_csv = Path(config.correlation_result_filename).resolve()
    output_dir = final_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = sanitize_prefix(config.tmp_name)

    # final_only means exactly one persistent pipeline output: final_csv.
    effective_save_before_bh = (
        config.save_before_bh_csv and not config.final_only
    )
    before_bh_csv = (
        derive_before_bh_path(final_csv)
        if effective_save_before_bh
        else None
    )
    _check_csv_collisions(
        before_bh_csv=before_bh_csv,
        final_csv=final_csv,
        overwrite=config.overwrite,
    )

    work_dir = _prepare_work_dir(
        output_dir=output_dir,
        prefix=prefix,
        final_only=config.final_only,
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

    merged_result = None
    bh_result = None
    completed = False
    try:
        print("Calculating and merging filtered/unfiltered correlations...")
        merged_result = calculate_merged_correlations(
            filtered_values,
            values,
            jobs=config.n_jobs,
            backend="cpu",
            block_size=config.block_size,
            output_dtype=config.output_dtype,
            output_dir=work_dir,
            prefix="merged",
            blas_threads=config.blas_threads,
            overwrite=False,
        )
        del filtered_values
        del values
        gc.collect()

        merged_corr = merged_result["corr"]
        merged_p = merged_result["p_value"]

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
        bh_metadata_path = work_dir / "bh_metadata.json"
        bh_result = adjust_pvalues_bh(
            merged_p,
            alpha=0.05,
            scan_chunk_size=config.bh_chunk_size,
            metadata_path=bh_metadata_path,
            show_progress=config.show_progress,
        )

        print(f"Writing adjusted-p < 0.05 CSV: {final_csv}")
        significant_pairs = write_bh_significant_csv(
            substance_ids=substance_ids,
            correlations=merged_corr,
            p_values=merged_p,
            sorted_significant_pvalues=(
                bh_result["sorted_significant_pvalues"]
            ),
            adjusted_significant=bh_result["adjusted_significant"],
            raw_cutoff=bh_result["raw_cutoff"],
            rejection_count=bh_result["rejection_count"],
            adjusted_threshold=0.05,
            output_path=final_csv,
            pair_chunk_size=config.pair_chunk_size,
            overwrite=config.overwrite,
            show_progress=config.show_progress,
        )
        completed = True

        retained_work_dir = None if config.final_only else work_dir
        retained_merged_metadata = (
            None
            if config.final_only
            else Path(merged_result["metadata_path"])
        )
        retained_bh_metadata = (
            None
            if config.final_only
            else Path(bh_result["metadata_path"])
        )
        return PipelineResult(
            before_bh_csv=before_bh_csv,
            significant_csv=final_csv,
            work_dir=retained_work_dir,
            merged_metadata=retained_merged_metadata,
            bh_metadata=retained_bh_metadata,
            total_pairs=total_pairs,
            significant_pairs=significant_pairs,
        )
    finally:
        # Release potentially tens of GiB before deleting temporary files.
        if bh_result is not None:
            bh_result.pop("sorted_significant_pvalues", None)
            bh_result.pop("adjusted_significant", None)
        if merged_result is not None:
            close_memmaps(
                merged_result.get("corr"),
                merged_result.get("p_value"),
            )
        gc.collect()

        # On successful final_only runs, leave only the requested final CSV.
        # On failure, retain the work directory for diagnosis/restart analysis.
        if config.final_only and completed and work_dir.exists():
            shutil.rmtree(work_dir)


def _validate_config(config: PipelineConfig) -> None:
    if config.compounds_num <= 0:
        raise ValueError("compounds_num must be positive")
    if config.samples_num <= 0:
        raise ValueError("samples_num must be positive")
    if config.n_jobs <= 0:
        raise ValueError("n_jobs must be positive")
    if config.pair_chunk_size <= 0:
        raise ValueError("pair_chunk_size must be positive")
    if config.bh_chunk_size <= 0:
        raise ValueError("bh_chunk_size must be positive")
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


def _prepare_work_dir(
    *,
    output_dir: Path,
    prefix: str,
    final_only: bool,
    overwrite: bool,
) -> Path:
    if final_only:
        work_dir = output_dir / f".{prefix}_work_{uuid.uuid4().hex[:10]}"
        work_dir.mkdir(parents=True, exist_ok=False)
        return work_dir

    work_dir = output_dir / f"{prefix}_work"
    if work_dir.exists():
        if not overwrite:
            raise FileExistsError(f"work directory already exists: {work_dir}")
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=False)
    return work_dir
