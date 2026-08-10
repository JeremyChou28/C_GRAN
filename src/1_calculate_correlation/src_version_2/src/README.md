# Pairwise-complete correlation pipeline

This project refactors the original `calculate_correlation.py` workflow into:

- `entry.py`: command-line entry only
- `engine.py`: blocked pairwise-complete Pearson calculation
- `storage.py`: packed-upper memmap and pair-order helpers
- `stats_utils.py`: legacy filtered/unfiltered merge logic and BH correction
- `io_utils.py`: TSV loading, NumPy IQR filtering, chunked CSV writing
- `pipeline.py`: end-to-end orchestration
- `tests/reference_legacy.py`: the original pandas/row-pair algorithm used as the test oracle

## Outputs

For `--correlation_result_filename output/result.csv`, the pipeline writes:

1. `output/result_before_bh.csv` (optional; enabled by default)
   - all non-diagonal pairs after merging filtered and unfiltered results
   - columns: `Substance 1`, `Substance 2`, `Correlation`, `P-Value`
   - pass `--skip_before_bh_csv` to avoid writing this very large file

2. `output/result.csv`
   - only pairs with `adjusted p-value < 0.05`
   - columns: `Substance 1`, `Substance 2`, `Correlation`, `P-Value`, `adjusted p-value`

Filtered, unfiltered, merged and adjusted-p memmaps are retained in the same output directory.

The engine returns only:

- `corr`
- `p_value`
- `common_count`
- metadata

There is no `illegal_mask`; `common_count < 3` is the corresponding condition.

## Installation

```bash
pip install -r requirements.txt
```

## Run

```bash
python entry.py \
  --intensity_file input/AREA_MSV000094018.txt \
  --compounds_num 36000 \
  --samples_num 30 \
  --correlation_result_filename output/correlation_results_MSV000094018_pos.csv \
  --n_jobs 8 \
  --tmp_name coral/MSV000094018_pos \
  --block_size 1024 \
  --pair_chunk_size 1000000
```

`tmp_name` is sanitized and used as the retained memmap filename prefix. For example, `coral/MSV000094018_pos` becomes `coral_MSV000094018_pos`.

Use `--overwrite` to replace an existing run with the same output paths and prefix.

The pre-BH CSV is enabled by default for backward compatibility. To skip it
while still running BH from the merged p-value memmap:

```bash
python entry.py \
  ... \
  --skip_before_bh_csv
```

You may explicitly select the default behavior with `--save_before_bh_csv`.
The two flags are mutually exclusive.

## Storage for 36,000 rows

There are 647,982,000 non-diagonal pairs and 648,018,000 packed-upper entries including the diagonal.

With float64 correlation and p-value storage:

- when `samples_num <= 255`, `common_count` uses uint8
- two engine result sets plus merged correlation, raw p-value and adjusted p-value use about **35 GiB** of memmap space
- when `common_count` requires uint16, the memmaps use about **36.2 GiB**

The two CSV files require additional disk space. The pre-BH CSV contains every non-diagonal pair and can be very large.

BH correction is intentionally performed as one global in-memory `statsmodels.stats.multitest.multipletests(..., method="fdr_bh")` call.

## Tests

```bash
pytest -q
```

The tests include:

- SciPy pairwise Pearson comparison
- `common_count` comparison
- packed/full equivalence
- memmap reopen and overwrite behavior
- NumPy IQR filtering versus the original pandas implementation
- legacy Fisher/NaN/min-|r| behavior comparison
- BH comparison against statsmodels
- complete old-versus-new pipeline comparison for both CSV outputs
- direct CLI execution
