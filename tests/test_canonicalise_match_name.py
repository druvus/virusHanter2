"""Unit tests for the BLAST-title canonicalisation helper."""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.functions import (  # noqa: E402
    _is_strain_like_name,
    canonicalise_blast_match_name,
)


# Minimal taxdump fragment matching the real NCBI structure for the
# EBV-1 / EBV-2 case the canonicalisation is designed to collapse.
# Lines follow the standard ``id\t|\tparent\t|\trank\t|\t...`` shape.
_NODES_DMP = textwrap.dedent(
    """\
    1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    10375\t|\t1\t|\tgenus\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    10376\t|\t10375\t|\tno rank\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    12509\t|\t10376\t|\tno rank\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    99999\t|\t1\t|\tno rank\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    """
)


# Each scientific-name row carries the four-column shape NCBI uses;
# only ``scientific name`` rows are consumed by the canonicaliser so
# unrelated common-name / acronym rows are present here too as a
# stress test of the row-type filter.
_NAMES_DMP = textwrap.dedent(
    """\
    10375\t|\tLymphocryptovirus\t|\t\t|\tscientific name\t|
    10376\t|\thuman gammaherpesvirus 4\t|\t\t|\tscientific name\t|
    10376\t|\tEBV\t|\t\t|\tacronym\t|
    10376\t|\tHuman herpesvirus 4\t|\t\t|\tequivalent name\t|
    12509\t|\tHuman herpesvirus 4 type 2\t|\t\t|\tscientific name\t|
    99999\t|\tSome species type 1\t|\t\t|\tscientific name\t|
    """
)


@pytest.fixture
def taxdump(tmp_path: Path) -> tuple[Path, Path]:
    nodes = tmp_path / "nodes.dmp"
    names = tmp_path / "names.dmp"
    nodes.write_text(_NODES_DMP)
    names.write_text(_NAMES_DMP)
    return nodes, names


@pytest.fixture
def parquet_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "name": [
                "NC_007605.1 Human gammaherpesvirus 4, complete genome",
                "NC_009334.1 Human herpesvirus 4, complete genome",
                "NC_999999.1 Some species type 1",
            ],
            "sequence": ["A", "C", "G"],
            "tax_id": [10376, 12509, 99999],
        }
    )


def test_is_strain_like_name_matches_only_with_word_boundaries():
    assert _is_strain_like_name("Human herpesvirus 4 type 2")
    assert _is_strain_like_name("HSV-1 strain F")
    assert _is_strain_like_name("Influenza A virus isolate X")
    # No marker -> not strain-like.
    assert not _is_strain_like_name("human gammaherpesvirus 4")
    assert not _is_strain_like_name("Lymphocryptovirus")
    # "type" embedded inside a longer word should not match — the
    # marker is " type " (space-padded).
    assert not _is_strain_like_name("Phenotype X virus")
    assert not _is_strain_like_name("")


def test_canonicalise_collapses_ebv_type_2_onto_species(parquet_df, taxdump):
    nodes, names = taxdump
    blastn = pd.DataFrame(
        {
            "name": ["contig_A", "contig_B"],
            "match_name": [
                "Human gammaherpesvirus 4, complete genome",
                "Human herpesvirus 4, complete genome",
            ],
            "accession": ["NC_007605", "NC_009334"],
            "read_len": [3000, 2500],
        }
    )
    out = canonicalise_blast_match_name(blastn, parquet_df, str(nodes), str(names))
    # Both rows now share the canonical species name.
    assert list(out["match_name"]) == [
        "human gammaherpesvirus 4",
        "human gammaherpesvirus 4",
    ]
    # The raw BLAST title is preserved for audit.
    assert list(out["match_name_raw"]) == [
        "Human gammaherpesvirus 4, complete genome",
        "Human herpesvirus 4, complete genome",
    ]


def test_canonicalise_keeps_unmapped_rows_unchanged(parquet_df, taxdump):
    nodes, names = taxdump
    blastn = pd.DataFrame(
        {
            "name": ["contig_X"],
            "match_name": ["Mystery virus X, complete genome"],
            "accession": ["NC_NOT_IN_PARQUET"],
            "read_len": [1000],
        }
    )
    out = canonicalise_blast_match_name(blastn, parquet_df, str(nodes), str(names))
    assert list(out["match_name"]) == ["Mystery virus X, complete genome"]
    assert list(out["match_name_raw"]) == ["Mystery virus X, complete genome"]


def test_canonicalise_walks_chain_of_strain_like_names(parquet_df, taxdump):
    nodes, names = taxdump
    # A taxid whose name *is* strain-like and whose parent is the
    # taxonomic root. The walk should stop at the input and use its
    # own name because no non-strain-like ancestor exists below
    # taxid 1.
    blastn = pd.DataFrame(
        {
            "name": ["contig_Y"],
            "match_name": ["Some species type 1"],
            "accession": ["NC_999999"],
            "read_len": [800],
        }
    )
    out = canonicalise_blast_match_name(blastn, parquet_df, str(nodes), str(names))
    # Parent is root (taxid 1); the walk stops with the original name
    # untouched rather than substituting "root" or similar.
    assert list(out["match_name"]) == ["Some species type 1"]
    assert list(out["match_name_raw"]) == ["Some species type 1"]


def test_canonicalise_no_taxdump_degrades_gracefully(parquet_df):
    blastn = pd.DataFrame(
        {
            "name": ["contig_A"],
            "match_name": ["Human gammaherpesvirus 4, complete genome"],
            "accession": ["NC_007605"],
            "read_len": [3000],
        }
    )
    out = canonicalise_blast_match_name(blastn, parquet_df, None, None)
    assert list(out["match_name"]) == ["Human gammaherpesvirus 4, complete genome"]
    assert list(out["match_name_raw"]) == [
        "Human gammaherpesvirus 4, complete genome"
    ]


def test_canonicalise_empty_blast_passes_through():
    out = canonicalise_blast_match_name(pd.DataFrame(), pd.DataFrame(), None, None)
    assert out.empty
