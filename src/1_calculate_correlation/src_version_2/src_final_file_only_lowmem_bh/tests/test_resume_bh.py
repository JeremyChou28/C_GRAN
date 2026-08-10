import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import PipelineConfig, run_pipeline


def test_resume_bh_cli_reuses_retained_merged_files(tmp_path):
    project = Path(__file__).resolve().parents[1]
    values = np.array(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            [2.0, 4.0, 6.0, 8.0, 10.0, 12.0],
            [6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            [1.0, 2.5, 2.0, 5.0, 4.0, 7.0],
            [1.0, 2.0, 3.0, 4.0, 5.0, 100.0],
        ]
    )
    ids = [f"R_{i}" for i in range(values.shape[0])]
    input_path = tmp_path / "input.tsv"
    pd.DataFrame(values, index=ids).to_csv(input_path, sep="\t")

    original_final = tmp_path / "original.csv"
    result = run_pipeline(
        PipelineConfig(
            intensity_file=input_path,
            compounds_num=values.shape[0],
            samples_num=values.shape[1],
            correlation_result_filename=original_final,
            n_jobs=1,
            tmp_name="resume-source",
            block_size=2,
            pair_chunk_size=3,
            bh_chunk_size=4,
            save_before_bh_csv=False,
            show_progress=False,
        )
    )
    resumed_final = tmp_path / "resumed.csv"
    completed = subprocess.run(
        [
            sys.executable,
            str(project / "resume_bh.py"),
            "--intensity_file",
            str(input_path),
            "--merged_metadata",
            str(result.merged_metadata),
            "--correlation_result_filename",
            str(resumed_final),
            "--pair_chunk_size",
            "3",
            "--bh_chunk_size",
            "4",
            "--no_progress",
        ],
        cwd=project,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr

    original = pd.read_csv(original_final)
    resumed = pd.read_csv(resumed_final)
    pd.testing.assert_frame_equal(original, resumed)
