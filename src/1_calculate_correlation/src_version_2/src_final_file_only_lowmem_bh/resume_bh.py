#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import gc
import shutil
from pathlib import Path

import numpy as np

from io_utils import load_substance_ids, write_bh_significant_csv
from stats_utils import adjust_pvalues_bh
from storage import close_memmaps, read_json


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resume the low-memory BH and final-CSV stages from an existing "
            "fused merged_corr.dat / merged_p_value.dat work directory."
        )
    )
    parser.add_argument(
        "--intensity_file",
        required=True,
        type=Path,
        help="original TSV; only the first substance-ID column is read",
    )
    parser.add_argument(
        "--merged_metadata",
        required=True,
        type=Path,
        help="path to the retained merged_metadata.json",
    )
    parser.add_argument(
        "--correlation_result_filename",
        required=True,
        type=Path,
        help="final CSV containing adjusted-p < 0.05 pairs",
    )
    parser.add_argument(
        "--pair_chunk_size",
        default=1_000_000,
        type=int,
        help="number of pairs formatted per CSV chunk",
    )
    parser.add_argument(
        "--bh_chunk_size",
        default=10_000_000,
        type=int,
        help="number of raw p-values scanned per low-memory BH chunk",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace an existing final CSV",
    )
    parser.add_argument(
        "--delete_work_dir_on_success",
        action="store_true",
        help="delete the directory containing merged_metadata.json after success",
    )
    parser.add_argument(
        "--no_progress",
        action="store_true",
        help="disable tqdm progress bars",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    metadata_path = args.merged_metadata.resolve()
    metadata = read_json(metadata_path)
    if not metadata.get("complete", False):
        raise ValueError(f"merged result is incomplete: {metadata_path}")
    if metadata.get("kind") != "fused_filtered_unfiltered_legacy_merge":
        raise ValueError("metadata is not a fused merged-correlation result")

    n_rows = int(metadata["n_rows"])
    pair_total = int(metadata["pair_count"])
    expected_pairs = n_rows * (n_rows - 1) // 2
    if pair_total != expected_pairs:
        raise ValueError(
            f"metadata pair_count mismatch: {pair_total} != {expected_pairs}"
        )

    substance_ids = load_substance_ids(
        args.intensity_file,
        expected_rows=n_rows,
    )
    paths = metadata["paths"]
    dtype = np.dtype(metadata.get("merged_output_dtype", "float64"))
    merged_corr = np.memmap(
        paths["corr"],
        dtype=dtype,
        mode="r",
        shape=(pair_total,),
    )
    merged_p = np.memmap(
        paths["p_value"],
        dtype=dtype,
        mode="r",
        shape=(pair_total,),
    )

    work_dir = metadata_path.parent
    final_csv = args.correlation_result_filename.resolve()
    if args.delete_work_dir_on_success:
        try:
            final_csv.relative_to(work_dir.resolve())
        except ValueError:
            pass
        else:
            raise ValueError(
                "final CSV must not be inside the work directory when "
                "--delete_work_dir_on_success is used"
            )

    bh_result = None
    completed = False
    try:
        bh_result = adjust_pvalues_bh(
            merged_p,
            alpha=0.05,
            scan_chunk_size=args.bh_chunk_size,
            metadata_path=work_dir / "bh_metadata.json",
            show_progress=not args.no_progress,
        )
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
            pair_chunk_size=args.pair_chunk_size,
            overwrite=args.overwrite,
            show_progress=not args.no_progress,
        )
        completed = True
        print(f"Final CSV: {final_csv}")
        print(f"Valid p-values: {bh_result['valid_count']:,}")
        print(f"BH rejection count (adjusted <= 0.05): {bh_result['rejection_count']:,}")
        print(f"Written pairs (adjusted < 0.05): {significant_pairs:,}")
        return 0
    finally:
        if bh_result is not None:
            bh_result.pop("sorted_significant_pvalues", None)
            bh_result.pop("adjusted_significant", None)
        close_memmaps(merged_corr, merged_p)
        gc.collect()
        if completed and args.delete_work_dir_on_success:
            shutil.rmtree(work_dir)


if __name__ == "__main__":
    raise SystemExit(main())
