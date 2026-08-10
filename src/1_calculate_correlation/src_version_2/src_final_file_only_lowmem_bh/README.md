# Pairwise-complete correlation pipeline — fused merge + low-memory exact BH

This project computes pairwise-complete Pearson correlations between rows,
keeps the original filtered/unfiltered merge semantics, applies one exact
global Benjamini-Hochberg correction, and writes the significant pair CSV.

The implementation is designed for very large row counts:

- filtered and unfiltered correlations are calculated and merged inside each
  in-memory block;
- only `merged_corr.dat` and `merged_p_value.dat` are used as large temporary
  files;
- the BH implementation does **not** sort every finite p-value;
- in `--final_only` mode, successful runs delete all work files and retain only
  the final CSV.

## Low-memory exact BH

The previous version called:

```python
multipletests(all_valid_pvalues, method="fdr_bh")
```

For billions of tests, statsmodels creates several arrays as large as the full
valid-p-value vector: a value copy, sort indices, sorted values, ranks,
adjusted values, and order-restoration arrays. This can require several hundred
GiB of RAM.

The new implementation computes exactly the part needed for the final output
`adjusted p-value < 0.05`:

1. sequentially count finite p-values (`m`);
2. find the exact BH rejection count with the monotone fixed point
   `R = count(p <= 0.05 * R / m)`;
3. extract and sort only the `R` rejected p-values;
4. calculate their exact BH-adjusted p-values;
5. scan the pair arrays once more and write only adjusted-p `< 0.05` rows.

This is not blockwise or approximate BH. Tests compare the selected hypotheses
and adjusted values directly against `statsmodels.stats.multitest.multipletests`.

Peak BH memory is approximately:

```text
16 * R bytes + BH chunk buffers
```

where `R` is the BH rejection count, rather than the total number of valid
p-values. The tradeoff is several sequential scans of `merged_p_value.dat`.

## Files

- `entry.py`: normal command-line entry
- `resume_bh.py`: resume BH/final CSV from retained merged files
- `pipeline.py`: end-to-end orchestration and work cleanup
- `engine.py`: blocked Pearson calculation and fused filtered/unfiltered merge
- `stats_utils.py`: legacy merge semantics and exact low-memory BH
- `storage.py`: memmap, packed-triangle and pair-order helpers
- `io_utils.py`: TSV loading, NumPy IQR filtering and chunked CSV output
- `tests/reference_legacy.py`: original algorithm used as the test oracle

## Installation

```bash
pip install -r requirements.txt
```

## Recommended run

```bash
python entry.py \
  --intensity_file ../input/coral_processed_area/AREA_MSV000082083.txt \
  --compounds_num 131620 \
  --samples_num 321 \
  --correlation_result_filename ../output/MSV000082083.csv \
  --n_jobs 32 \
  --tmp_name MSV000082083 \
  --block_size 512 \
  --pair_chunk_size 1000000 \
  --bh_chunk_size 10000000 \
  --final_only
```

`--bh_chunk_size` controls only sequential BH scans. A value of ten million
uses roughly tens to a few hundred MiB of temporary memory per scan and avoids
excessive Python-loop overhead.

With `--final_only`:

- the complete pre-BH CSV is skipped;
- a hidden temporary directory is created beside the final CSV;
- only `merged_corr.dat` and `merged_p_value.dat` are used as large temporary
  files;
- after the final CSV is successfully written, the temporary directory is
  deleted;
- if the process fails, the directory is retained for diagnosis and resume.

## Resume the run that failed during BH

A BH-stage failure occurs after correlation and fused merge have completed. Do
not delete the retained work directory. Locate its metadata:

```bash
find ../output -name merged_metadata.json -print
```

Then run:

```bash
python resume_bh.py \
  --intensity_file ../input/coral_processed_area/AREA_MSV000082083.txt \
  --merged_metadata ../output/.MSV000082083_work_xxxxxxxxxx/merged_metadata.json \
  --correlation_result_filename ../output/MSV000082083.csv \
  --pair_chunk_size 1000000 \
  --bh_chunk_size 10000000
```

This reads only the substance-ID column from the original TSV and reuses the
existing `merged_corr.dat` and `merged_p_value.dat`. It does not recalculate
correlations.

After verifying the output, the work directory can be removed manually. Or add:

```bash
--delete_work_dir_on_success
```

The final CSV must be outside the work directory when that option is used.

## Retain work files during a new run

Omit `--final_only` to retain:

```text
output/<sanitized_tmp_name>_work/
├── merged_corr.dat
├── merged_p_value.dat
├── merged_metadata.json
└── bh_metadata.json
```

The full pre-BH CSV is enabled by default in this mode. Skip it with:

```bash
--skip_before_bh_csv
```

The final CSV always contains only pairs satisfying:

```text
adjusted p-value < 0.05
```

with columns:

```text
Substance 1, Substance 2, Correlation, P-Value, adjusted p-value
```

## Disk use for 131,620 rows

The non-diagonal pair count is:

```text
8,661,846,390
```

The two float64 merged arrays require approximately:

- `merged_corr.dat`: 64.54 GiB
- `merged_p_value.dat`: 64.54 GiB
- total: 129.07 GiB

The final CSV requires additional space. In `--final_only` mode these temporary
files are deleted after successful CSV generation.

## Tests

```bash
pytest -q
```

The suite covers:

- SciPy pairwise Pearson comparison;
- common-count validation for the general engine;
- packed/full equivalence;
- NumPy IQR filtering versus the old pandas implementation;
- legacy Fisher, NaN and min-|r| behavior;
- fused block output versus the original algorithm;
- exact low-memory BH versus statsmodels, including ties and NaNs;
- both CSV outputs versus the original pipeline;
- `--final_only` cleanup behavior;
- direct `python entry.py ...` execution;
- resuming BH from retained merged files.
