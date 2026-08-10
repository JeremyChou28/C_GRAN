#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from pipeline import PipelineConfig, run_pipeline


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate pairwise-complete row correlations, merge filtered and "
            "unfiltered results with the legacy logic, then apply global BH."
        )
    )
    parser.add_argument(
        "--intensity_file",
        required=True,
        type=Path,
        help="input TSV file; first column is the substance ID",
    )
    parser.add_argument(
        "--compounds_num",
        required=True,
        type=int,
        help="expected number of compounds/rows",
    )
    parser.add_argument(
        "--samples_num",
        required=True,
        type=int,
        help="expected number of samples/features",
    )
    parser.add_argument(
        "--correlation_result_filename",
        default=Path("correlation_results.csv"),
        type=Path,
        help=(
            "final CSV containing pairs with adjusted p-value < 0.05"
        ),
    )
    parser.add_argument(
        "--n_jobs",
        default=32,
        type=int,
        help="number of block workers",
    )
    parser.add_argument(
        "--tmp_name",
        required=True,
        type=str,
        help="prefix for the work directory and retained work files",
    )
    parser.add_argument(
        "--block_size",
        default=None,
        type=int,
        help="row block size; selected automatically when omitted",
    )
    parser.add_argument(
        "--output_dtype",
        choices=("float32", "float64"),
        default="float64",
        help="storage dtype for engine corr and p-value memmaps",
    )
    parser.add_argument(
        "--blas_threads",
        default=None,
        type=int,
        help="BLAS threads available to each numerical call",
    )
    parser.add_argument(
        "--pair_chunk_size",
        default=1_000_000,
        type=int,
        help="number of pairs processed per CSV chunk",
    )
    parser.add_argument(
        "--bh_chunk_size",
        default=10_000_000,
        type=int,
        help=(
            "number of raw p-values scanned per low-memory BH chunk; "
            "larger values reduce loop overhead but use more temporary RAM"
        ),
    )
    before_bh_group = parser.add_mutually_exclusive_group()
    before_bh_group.add_argument(
        "--save_before_bh_csv",
        dest="save_before_bh_csv",
        action="store_true",
        help=(
            "write the complete merged CSV before BH correction "
            "(default behavior)"
        ),
    )
    before_bh_group.add_argument(
        "--skip_before_bh_csv",
        dest="save_before_bh_csv",
        action="store_false",
        help=(
            "skip the complete pre-BH CSV; BH still runs directly from "
            "the merged p-value memmap"
        ),
    )
    parser.set_defaults(save_before_bh_csv=True)
    parser.add_argument(
        "--final_only",
        action="store_true",
        help=(
            "persist only the final adjusted-p < 0.05 CSV; this also skips "
            "the pre-BH CSV and deletes the temporary merged memmaps after "
            "a successful run"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace existing CSV, memmap and metadata files",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    start = time.time()

    result = run_pipeline(
        PipelineConfig(
            intensity_file=args.intensity_file,
            compounds_num=args.compounds_num,
            samples_num=args.samples_num,
            correlation_result_filename=args.correlation_result_filename,
            n_jobs=args.n_jobs,
            tmp_name=args.tmp_name,
            block_size=args.block_size,
            output_dtype=np.dtype(args.output_dtype),
            blas_threads=args.blas_threads,
            pair_chunk_size=args.pair_chunk_size,
            bh_chunk_size=args.bh_chunk_size,
            save_before_bh_csv=args.save_before_bh_csv,
            final_only=args.final_only,
            overwrite=args.overwrite,
        )
    )

    print("Pipeline completed.")
    if result.before_bh_csv is None:
        print("Pre-BH CSV: skipped")
    else:
        print(f"Pre-BH CSV: {result.before_bh_csv}")
    print(f"Adjusted-p < 0.05 CSV: {result.significant_csv}")
    if result.work_dir is None:
        print("Work files: removed")
    else:
        print(f"Work files: {result.work_dir}")
    print(f"Total non-diagonal pairs: {result.total_pairs:,}")
    print(f"Significant pairs: {result.significant_pairs:,}")
    print(f"Spend time: {time.time() - start:.3f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
