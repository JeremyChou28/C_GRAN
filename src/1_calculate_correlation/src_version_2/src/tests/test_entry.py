import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def test_entry_runs_as_direct_script(tmp_path):
    project = Path(__file__).resolve().parents[1]
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
    final_path = tmp_path / "final.csv"

    command = [
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
    completed = subprocess.run(
        command,
        cwd=project,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    assert final_path.exists()
    assert (tmp_path / "final_before_bh.csv").exists()
    assert "Pipeline completed." in completed.stdout
    assert (tmp_path / "entry_test_filtered_metadata.json").exists()
    assert (tmp_path / "entry_test_unfiltered_common_count.dat").exists()



def test_entry_can_skip_before_bh_csv(tmp_path):
    project = Path(__file__).resolve().parents[1]
    values = np.array(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [2.0, 4.0, 6.0, 8.0, 10.0],
            [5.0, 4.0, 3.0, 2.0, 1.0],
            [1.0, 3.0, 2.0, 5.0, 4.0],
        ]
    )
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=["A", "B", "C", "D"]).to_csv(
        input_path,
        sep="\t",
    )
    final_path = tmp_path / "final.csv"

    command = [
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
        "entry-skip",
        "--block_size",
        "2",
        "--pair_chunk_size",
        "2",
        "--skip_before_bh_csv",
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
