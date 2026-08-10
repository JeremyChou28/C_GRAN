import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _write_input(tmp_path):
    values = np.array(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [2.0, 4.0, 6.0, 8.0, 10.0],
            [5.0, 4.0, 3.0, 2.0, 1.0],
            [1.0, 3.0, 2.0, 5.0, 4.0],
        ]
    )
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=["A_1", "B_2", "C_3", "D_4"]).to_csv(
        input_path,
        sep="\t",
    )
    return input_path


def _base_command(project, input_path, final_path):
    return [
        sys.executable,
        str(project / "entry.py"),
        "--intensity_file",
        str(input_path),
        "--compounds_num",
        "4",
        "--samples_num",
        "5",
        "--correlation_result_filename",
        str(final_path),
        "--n_jobs",
        "1",
        "--tmp_name",
        "entry/test",
        "--block_size",
        "2",
        "--pair_chunk_size",
        "2",
    ]


def test_entry_runs_as_direct_script_and_retains_fused_work(tmp_path):
    project = Path(__file__).resolve().parents[1]
    input_path = _write_input(tmp_path)
    final_path = tmp_path / "final.csv"

    completed = subprocess.run(
        _base_command(project, input_path, final_path),
        cwd=project,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    assert final_path.exists()
    assert (tmp_path / "final_before_bh.csv").exists()
    assert "Pipeline completed." in completed.stdout
    work_dir = tmp_path / "entry_test_work"
    assert work_dir.exists()
    assert (work_dir / "merged_corr.dat").exists()
    assert (work_dir / "merged_p_value.dat").exists()
    assert not any("filtered" in p.name for p in work_dir.iterdir())


def test_entry_final_only_leaves_only_final_csv(tmp_path):
    project = Path(__file__).resolve().parents[1]
    input_path = _write_input(tmp_path)
    final_path = tmp_path / "final.csv"
    command = _base_command(project, input_path, final_path) + ["--final_only"]

    completed = subprocess.run(
        command,
        cwd=project,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    assert final_path.exists()
    assert not (tmp_path / "final_before_bh.csv").exists()
    assert "Pre-BH CSV: skipped" in completed.stdout
    assert "Work files: removed" in completed.stdout
    assert {p.name for p in tmp_path.iterdir()} == {"input.tsv", "final.csv"}


def test_entry_can_skip_before_bh_csv_but_retain_work(tmp_path):
    project = Path(__file__).resolve().parents[1]
    input_path = _write_input(tmp_path)
    final_path = tmp_path / "final.csv"
    command = _base_command(project, input_path, final_path) + [
        "--skip_before_bh_csv"
    ]

    completed = subprocess.run(
        command,
        cwd=project,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    assert final_path.exists()
    assert not (tmp_path / "final_before_bh.csv").exists()
    assert "Pre-BH CSV: skipped" in completed.stdout
    assert (tmp_path / "entry_test_work").exists()
