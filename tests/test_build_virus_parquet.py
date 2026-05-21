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

from scripts.build_virus_parquet import (  # noqa: E402
    enrich_with_taxdump,
    find_genus_taxid,
    parse_nodes_dmp,
    pick_longest_per_taxid,
)


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


# ---------------------------------------------------------------------------
# parse_nodes_dmp / find_genus_taxid / enrich_with_taxdump
# ---------------------------------------------------------------------------


def _write_nodes_dmp(tmp_path, rows: list[tuple[int, int, str]]) -> Path:
    """Write a tiny synthetic nodes.dmp with NCBI's ``\\t|\\t`` delimiter."""
    path = tmp_path / "nodes.dmp"
    with open(path, "w") as fh:
        for tid, parent, rank in rows:
            fh.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\n")
    return path


def test_parse_nodes_dmp_basic(tmp_path):
    path = _write_nodes_dmp(
        tmp_path,
        [
            (1, 1, "no rank"),
            (10239, 1, "superkingdom"),
            (687331, 10239, "genus"),
            (3048424, 687331, "species"),
        ],
    )
    nodes = parse_nodes_dmp(path)
    assert nodes[1] == (1, "no rank")
    assert nodes[10239] == (1, "superkingdom")
    assert nodes[687331] == (10239, "genus")
    assert nodes[3048424] == (687331, "species")


def test_parse_nodes_dmp_missing_file(tmp_path):
    assert parse_nodes_dmp(tmp_path / "does-not-exist.dmp") == {}


def test_find_genus_taxid_walks_to_genus():
    nodes = {
        3048424: (687331, "species"),
        687331: (10239, "genus"),
        10239: (1, "superkingdom"),
        1: (1, "no rank"),
    }
    assert find_genus_taxid(3048424, nodes) == 687331


def test_find_genus_taxid_returns_zero_when_no_genus():
    nodes = {
        100: (10, "family"),
        10: (1, "order"),
        1: (1, "no rank"),
    }
    assert find_genus_taxid(100, nodes) == 0


def test_find_genus_taxid_respects_depth_limit():
    # Build a long synthetic chain (50 levels) with no genus. Depth
    # limit (default 20) must prevent unbounded walks.
    nodes = {i: (i - 1 if i > 1 else 1, "no rank") for i in range(1, 50)}
    assert find_genus_taxid(49, nodes) == 0


def test_find_genus_taxid_handles_cycle():
    nodes = {10: (20, "species"), 20: (10, "species")}
    assert find_genus_taxid(10, nodes) == 0


def test_enrich_with_taxdump_adds_columns(tmp_path):
    nodes = {
        3048424: (687331, "species"),
        687331: (10239, "genus"),
        10376: (10239, "species"),
        10239: (1, "superkingdom"),
        1: (1, "no rank"),
    }
    df = pd.DataFrame(
        {
            "name": ["a", "b", "c"],
            "sequence": ["AAAA", "TTTT", "GGGG"],
            "tax_id": [3048424, 10376, 999999],
        }
    )
    out = enrich_with_taxdump(df, nodes)
    assert list(out.columns) == [
        "name",
        "sequence",
        "tax_id",
        "rank",
        "genus_taxid",
    ]
    assert out.loc[out["tax_id"] == 3048424, "rank"].iloc[0] == "species"
    assert int(out.loc[out["tax_id"] == 3048424, "genus_taxid"].iloc[0]) == 687331
    # 10376 (a species under a superkingdom — no intermediate genus
    # in the synthetic lineage) walks up and finds no genus.
    assert int(out.loc[out["tax_id"] == 10376, "genus_taxid"].iloc[0]) == 0
    # 999999 is absent from the nodes dict.
    assert out.loc[out["tax_id"] == 999999, "rank"].iloc[0] == "unknown"
    assert int(out.loc[out["tax_id"] == 999999, "genus_taxid"].iloc[0]) == 0


def test_enrich_genus_self_when_row_already_at_genus_rank():
    """A row whose own tax_id is already at rank `genus` records
    itself as the genus_taxid rather than walking further up.
    """
    nodes = {
        687331: (10239, "genus"),
        10239: (1, "superkingdom"),
        1: (1, "no rank"),
    }
    df = pd.DataFrame(
        {"name": ["x"], "sequence": ["AA"], "tax_id": [687331]}
    )
    out = enrich_with_taxdump(df, nodes)
    assert int(out["genus_taxid"].iloc[0]) == 687331
