"""Concatenate per-sample `<sample>.per_virus.csv` files into a single
batch-level `per_virus_<batch>.csv`.

Pure concatenation: every row already carries `run_name`, `sample_name`,
and `date`, so no joining is needed. Empty inputs (samples whose Kraken
report had zero viral hits) contribute zero rows; the resulting batch
file still always has the documented columns.

Usage:

    python scripts/aggregate_per_virus.py \\
        --in <a.per_virus.csv> <b.per_virus.csv> ... \\
        --out per_virus_<batch>.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PER_VIRUS_COLUMNS = [
    "run_name",
    "sample_name",
    "date",
    "virus_name_kraken",
    "virus_taxid",
    "virus_name_kaiju",
    "contigs",
    "virus_reads_kraken2",
    "other_reads",
    "total_reads",
    "human_reads",
    "human_reads_percent",
    "non_human_reads",
    "non_human_reads_percent",
    "note",
    "specific_virus_rpm",
    "all_virus_rpm",
    "Completeness (% >5X)",
    "bases_above_5x",
    "mean_coverage",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inputs", nargs="+", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    frames: list[pd.DataFrame] = []
    for path in args.inputs:
        if not path.exists() or path.stat().st_size == 0:
            continue
        try:
            frames.append(pd.read_csv(path))
        except pd.errors.EmptyDataError:
            # A sample whose Kraken report had zero viral hits writes an
            # empty frame with no header row, which is a non-zero-byte
            # file the size guard above does not catch. Treat it like an
            # empty input; the canonical schema is restored below.
            continue
    if frames:
        df = pd.concat(frames, ignore_index=True)
        # Reindex to the canonical column order so the batch file has a
        # stable schema regardless of which sample wrote first. Any
        # additive trailing columns the per-sample step emits (per-
        # assembler `*_contigs`, geNomad scores) are preserved after the
        # fixed schema in their existing order, rather than dropped.
        extra = [c for c in df.columns if c not in PER_VIRUS_COLUMNS]
        df = df.reindex(columns=PER_VIRUS_COLUMNS + extra)
    else:
        df = pd.DataFrame(columns=PER_VIRUS_COLUMNS)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
