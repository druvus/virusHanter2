"""Unit tests for scripts/compare_parquet_kraken2."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.compare_parquet_kraken2 import (  # noqa: E402
    build_overlap_frame,
    parse_kraken_inspect,
    summarise,
)


def test_parse_kraken_inspect_extracts_taxonomy_ids():
    sample = (
        "# kraken2-inspect header line is dropped\n"
        "1.00\t100\t90\tD\t10239\t  Viruses\n"
        "0.50\t50\t40\tF\t11906\t    Retroviridae\n"
        "0.25\t25\t20\tS\t11908\t      Human T-lymphotropic virus\n"
        "bogus line\n"  # filtered out (< 6 fields)
        "\n"
    )
    out = parse_kraken_inspect(sample)
    assert out == {
        10239: "Viruses",
        11906: "Retroviridae",
        11908: "Human T-lymphotropic virus",
    }


def test_build_overlap_frame_classifies_into_three_buckets():
    parquet = {10, 20, 30}
    kraken = {20: "B", 30: "C", 40: "D"}
    out = build_overlap_frame(parquet, kraken)
    rows = out.set_index("tax_id").to_dict(orient="index")

    assert rows[10] == {
        "name": "",
        "in_parquet": True,
        "in_kraken2": False,
        "only_in": "parquet",
    }
    assert rows[20] == {
        "name": "B",
        "in_parquet": True,
        "in_kraken2": True,
        "only_in": "both",
    }
    assert rows[30] == {
        "name": "C",
        "in_parquet": True,
        "in_kraken2": True,
        "only_in": "both",
    }
    assert rows[40] == {
        "name": "D",
        "in_parquet": False,
        "in_kraken2": True,
        "only_in": "kraken2",
    }


def test_summarise_counts_match_set_arithmetic(tmp_path):
    parquet = {1, 2, 3, 4}
    kraken = {3: "x", 4: "y", 5: "z"}
    summary = summarise(parquet, kraken, tmp_path / "db")
    assert summary["parquet_taxids"] == 4
    assert summary["kraken2_db_taxids"] == 3
    assert summary["intersection_count"] == 2  # {3, 4}
    assert summary["parquet_only_count"] == 2  # {1, 2}
    assert summary["kraken2_only_count"] == 1  # {5}
