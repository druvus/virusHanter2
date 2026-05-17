"""Combine multiple virusHanter2 result folders into master CSVs.

Each Illumina run produces its own result folder (one per `SAMPLES`
input). After several batches have been processed the collaborator
wants a single file across runs. This script does not call Snakemake;
it only reads per-batch CSVs and concatenates them.

Each `--result-folder` argument is one batch folder, e.g.
`/path/to/RESULTS_FOLDER/<batch_id>/`, containing at least
`run_information_<batch_id>.csv` and (after the per-virus
addition) `per_virus_<batch_id>.csv`. The script picks them up by
glob.

Usage:

    python scripts/merge_runs.py \\
        --result-folder /path/to/results/<batch1> \\
        --result-folder /path/to/results/<batch2> \\
        --out-dir /path/to/master/
    # writes:
    #   <out-dir>/master_per_sample.csv
    #   <out-dir>/master_per_virus.csv
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd


def _collect(result_folder: Path, prefix: str) -> list[Path]:
    return sorted(result_folder.glob(f"{prefix}_*.csv"))


def merge_csvs(paths: list[Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--result-folder",
        dest="result_folders",
        required=True,
        action="append",
        type=Path,
        help=(
            "Path to one batch result folder. Repeat for each run to "
            "include."
        ),
    )
    p.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Directory to write master_per_sample.csv and master_per_virus.csv into.",
    )
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING"])
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(message)s")

    per_sample_paths: list[Path] = []
    per_virus_paths: list[Path] = []
    for folder in args.result_folders:
        folder = folder.resolve()
        if not folder.is_dir():
            logging.warning("Skipping non-directory: %s", folder)
            continue
        sample_csvs = _collect(folder, "run_information")
        virus_csvs = _collect(folder, "per_virus")
        logging.info(
            "%s: %d run_information_*.csv, %d per_virus_*.csv",
            folder,
            len(sample_csvs),
            len(virus_csvs),
        )
        per_sample_paths.extend(sample_csvs)
        per_virus_paths.extend(virus_csvs)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    sample_df = merge_csvs(per_sample_paths)
    virus_df = merge_csvs(per_virus_paths)

    sample_out = args.out_dir / "master_per_sample.csv"
    virus_out = args.out_dir / "master_per_virus.csv"
    sample_df.to_csv(sample_out, index=False)
    virus_df.to_csv(virus_out, index=False)

    logging.info("Wrote %s (%d rows)", sample_out, len(sample_df))
    logging.info("Wrote %s (%d rows)", virus_out, len(virus_df))


if __name__ == "__main__":
    main()
