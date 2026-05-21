"""Unit tests for scripts/build_virus_parquet.

Covers the helpers that drive the parquet rebuild: per-taxid
longest-sequence selection, schema invariance, and the
``tax_id == 0`` dropping when ``--one-rep-per-taxid`` is on.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_virus_parquet import pick_longest_per_taxid  # noqa: E402


def _frame(rows: list[tuple[str, str, int]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["name", "sequence", "tax_id"])


def test_pick_longest_per_taxid_keeps_one_row_per_taxid():
    df = _frame(
        [
            ("acc1 strain a", "AAAA", 100),
            ("acc2 strain b", "AAAAAAAA", 100),
            ("acc3 strain c", "AAA", 100),
            ("acc4 strain d", "AAAAA", 200),
            ("acc5 strain e", "AAAAAAA", 200),
        ]
    )
    out = pick_longest_per_taxid(df)
    assert out["tax_id"].tolist() == [100, 200]
    assert out.loc[out["tax_id"] == 100, "sequence"].iloc[0] == "AAAAAAAA"
    assert out.loc[out["tax_id"] == 200, "sequence"].iloc[0] == "AAAAAAA"


def test_pick_longest_per_taxid_drops_zero_taxid():
    df = _frame(
        [
            ("unmapped", "AAAAAAA", 0),
            ("mapped", "AAAA", 100),
        ]
    )
    out = pick_longest_per_taxid(df)
    assert (out["tax_id"] != 0).all()
    assert out["tax_id"].tolist() == [100]


def test_pick_longest_per_taxid_ties_resolve_by_input_order():
    df = _frame(
        [
            ("first", "AAAA", 100),
            ("second", "AAAA", 100),
            ("third", "AAAA", 100),
        ]
    )
    out = pick_longest_per_taxid(df)
    assert len(out) == 1
    assert out["name"].iloc[0] == "first"


def test_pick_longest_per_taxid_schema_invariant():
    df = _frame(
        [
            ("a", "AAAA", 100),
            ("b", "AA", 200),
        ]
    )
    out = pick_longest_per_taxid(df)
    assert list(out.columns) == ["name", "sequence", "tax_id"]


def test_pick_longest_per_taxid_empty_input():
    df = _frame([])
    out = pick_longest_per_taxid(df)
    assert out.empty
    assert list(out.columns) == ["name", "sequence", "tax_id"]


def test_pick_longest_per_taxid_only_zero_taxids():
    df = _frame([("a", "AA", 0), ("b", "AAA", 0)])
    out = pick_longest_per_taxid(df)
    assert out.empty
