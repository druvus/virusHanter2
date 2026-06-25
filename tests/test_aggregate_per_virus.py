"""Regression tests for the batch per-virus aggregation.

`scripts/aggregate_per_virus.py` concatenates the per-sample
`per_virus_*.csv` files into a single batch deliverable. The batch
schema must (a) carry the completeness column under the exact name the
per-sample step emits, and (b) preserve the additive trailing columns
(per-assembler contig counts, geNomad scores) rather than dropping them.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "aggregate_per_virus.py"

# The fixed leading columns the per-sample step emits, in order.
LEADING = [
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


def _sample_row(sample_name: str, completeness: float, **extra: object) -> dict[str, object]:
    row: dict[str, object] = {col: "" for col in LEADING}
    row["sample_name"] = sample_name
    row["Completeness (% >5X)"] = completeness
    row.update(extra)
    return row


def _run(inputs: list[Path], out: Path) -> pd.DataFrame:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--in", *map(str, inputs), "--out", str(out)],
        check=True,
    )
    return pd.read_csv(out)


def test_completeness_column_survives(tmp_path: Path) -> None:
    a = tmp_path / "per_virus_s1.csv"
    pd.DataFrame([_sample_row("s1", 87.5)]).to_csv(a, index=False)

    out = _run([a], tmp_path / "batch.csv")

    assert "Completeness (% >5X)" in out.columns
    assert "completeness_5x" not in out.columns
    # The real value must survive, not be reindexed away to NaN.
    assert out["Completeness (% >5X)"].iloc[0] == 87.5


def test_trailing_additive_columns_preserved(tmp_path: Path) -> None:
    a = tmp_path / "per_virus_s1.csv"
    pd.DataFrame(
        [
            _sample_row(
                "s1",
                50.0,
                megahit_contigs=3,
                metaspades_contigs=1,
                genomad_viral_contigs=2,
                genomad_max_virus_score=0.91,
            )
        ]
    ).to_csv(a, index=False)

    out = _run([a], tmp_path / "batch.csv")

    for col in ("megahit_contigs", "metaspades_contigs", "genomad_viral_contigs"):
        assert col in out.columns, f"{col} dropped by aggregation"
    assert out["megahit_contigs"].iloc[0] == 3
    assert out["genomad_max_virus_score"].iloc[0] == 0.91


def test_leading_schema_order_is_stable(tmp_path: Path) -> None:
    a = tmp_path / "per_virus_s1.csv"
    pd.DataFrame([_sample_row("s1", 10.0, megahit_contigs=1)]).to_csv(a, index=False)

    out = _run([a], tmp_path / "batch.csv")

    assert list(out.columns)[: len(LEADING)] == LEADING
