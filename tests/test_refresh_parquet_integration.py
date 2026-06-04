"""Integration tests for the refresh/refresh_virus_parquet.smk pipeline logic.

These tests exercise ``scripts/build_virus_parquet.py`` end-to-end with
synthetic FASTA and accession2taxid fixtures so the core parquet-build path
is covered without any NCBI network access.  A companion ``conftest.py``
provides the shared synthetic fixtures.

Run with::

    pytest tests/ -v

from the repository root.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_virus_parquet import build, enrich_with_taxdump, parse_nodes_dmp  # noqa: E402


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def virus_fasta(tmp_path: Path) -> Path:
    """Write a tiny synthetic viral nucleotide FASTA to a temp file."""
    fasta = tmp_path / "viral.fna"
    records = [
        ("NC_000001.1", "alpha virus isolate A", "ATCGATCGATCG" * 10),
        ("NC_000002.1", "beta virus strain B",   "GCTAGCTAGCTA" * 10),
        ("NC_000003.1", "gamma virus variant C",  "TTAATTAATTAA" * 10),
        # Second record for taxid 100001 — shorter than the first; the
        # one-rep-per-taxid filter should discard this one.
        ("NC_000004.1", "alpha virus isolate D",  "ATCG" * 3),
    ]
    with open(fasta, "w") as fh:
        for acc, desc, seq in records:
            fh.write(f">{acc} {desc}\n{seq}\n")
    return fasta


@pytest.fixture()
def accession2taxid_gz(tmp_path: Path) -> Path:
    """Write a minimal NCBI nucl_gb.accession2taxid.gz fixture.

    Format: accession TAB accession.version TAB taxid TAB gi
    """
    gz_path = tmp_path / "nucl_gb.accession2taxid.gz"
    rows = [
        ("NC_000001", "NC_000001.1", "100001", "1"),
        ("NC_000002", "NC_000002.1", "100002", "2"),
        ("NC_000003", "NC_000003.1", "100003", "3"),
        ("NC_000004", "NC_000004.1", "100001", "4"),  # duplicate taxid, shorter seq
    ]
    with gzip.open(gz_path, "wt") as fh:
        fh.write("accession\taccession.version\ttaxid\tgi\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    return gz_path


@pytest.fixture()
def nodes_dmp(tmp_path: Path) -> Path:
    """Write a minimal NCBI nodes.dmp covering the synthetic taxids."""
    path = tmp_path / "nodes.dmp"
    rows = [
        (1, 1, "no rank"),
        (10239, 1, "superkingdom"),
        (200000, 10239, "genus"),
        (100001, 200000, "species"),
        (100002, 200000, "species"),
        (100003, 200000, "species"),
    ]
    with open(path, "w") as fh:
        for tid, parent, rank in rows:
            fh.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\n")
    return path


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_build_parquet_basic_schema(virus_fasta: Path, accession2taxid_gz: Path) -> None:
    """build() returns a DataFrame with the expected three-column schema."""
    df, stats = build(
        fastas=[virus_fasta],
        taxid_path=accession2taxid_gz,
        one_rep_per_taxid=False,
    )
    assert set(df.columns) >= {"name", "sequence", "tax_id"}
    assert len(df) == 4
    assert (df["tax_id"] != 0).all()


def test_build_one_rep_per_taxid(virus_fasta: Path, accession2taxid_gz: Path) -> None:
    """One-rep-per-taxid mode keeps only the longest sequence per taxid."""
    df, stats = build(
        fastas=[virus_fasta],
        taxid_path=accession2taxid_gz,
        one_rep_per_taxid=True,
    )
    # Three distinct taxids; the shorter NC_000004 duplicate is dropped.
    assert len(df) == 3
    assert df["tax_id"].nunique() == 3
    # The longer NC_000001 record must survive over NC_000004.
    alpha_row = df.loc[df["tax_id"] == 100001]
    assert len(alpha_row) == 1
    assert "NC_000001" in alpha_row.iloc[0]["name"]


def test_build_stats_json_fields(virus_fasta: Path, accession2taxid_gz: Path) -> None:
    """build() returns a stats dict containing the expected headline keys."""
    _df, stats = build(
        fastas=[virus_fasta],
        taxid_path=accession2taxid_gz,
        one_rep_per_taxid=True,
    )
    assert "input_records" in stats
    assert "output_records" in stats
    assert "unique_taxids" in stats
    assert stats["output_records"] == stats["unique_taxids"] == 3


def test_build_stats_json_written_to_disk(
    tmp_path: Path, virus_fasta: Path, accession2taxid_gz: Path
) -> None:
    """build() + main() write build_stats.json next to the parquet."""
    import argparse
    from scripts.build_virus_parquet import main as _main
    import sys as _sys

    out_parquet = tmp_path / "out.parquet"
    _sys.argv = [
        "build_virus_parquet.py",
        "--fasta", str(virus_fasta),
        "--taxid", str(accession2taxid_gz),
        "--out", str(out_parquet),
        "--source", "refseq",
        "--one-rep-per-taxid",
        "--log-level", "WARNING",
    ]
    _main()

    assert out_parquet.exists(), "parquet file not written"
    stats_path = out_parquet.with_name("out_build_stats.json")
    assert stats_path.exists(), "build_stats.json sidecar not written"

    with open(stats_path) as fh:
        stats = json.load(fh)
    assert "build_date_utc" in stats
    assert stats["source"] == "refseq"
    assert stats["output_records"] == 3


def test_build_without_taxid_sets_zero(virus_fasta: Path) -> None:
    """When no accession2taxid is provided, every row gets tax_id == 0."""
    df, _stats = build(
        fastas=[virus_fasta],
        taxid_path=None,
        one_rep_per_taxid=False,
    )
    assert (df["tax_id"] == 0).all()


def test_build_with_taxdump_enrichment(
    virus_fasta: Path, accession2taxid_gz: Path, nodes_dmp: Path
) -> None:
    """When a nodes.dmp is supplied, rank and genus_taxid columns are added."""
    df, _stats = build(
        fastas=[virus_fasta],
        taxid_path=accession2taxid_gz,
        one_rep_per_taxid=True,
        taxdump_nodes_path=nodes_dmp,
    )
    assert "rank" in df.columns
    assert "genus_taxid" in df.columns
    species_rows = df.loc[df["rank"] == "species"]
    assert len(species_rows) == 3
    # All three species share genus 200000.
    assert (species_rows["genus_taxid"] == 200000).all()


def test_build_warns_on_missing_taxids(
    tmp_path: Path, accession2taxid_gz: Path
) -> None:
    """Records whose accession is absent from accession2taxid get tax_id 0."""
    fasta = tmp_path / "unmapped.fna"
    fasta.write_text(">UNMAPPED_ACC.1 synthetic unmapped\nATCGATCGATCG\n")
    df, stats = build(
        fastas=[fasta],
        taxid_path=accession2taxid_gz,
        one_rep_per_taxid=False,
    )
    assert len(df) == 1
    assert df.iloc[0]["tax_id"] == 0
    assert stats["unresolved_taxid_records"] == 1
