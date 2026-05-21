"""Unit tests for scripts/reformat_kaiju_headers."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.reformat_kaiju_headers import (  # noqa: E402
    collect_wanted_accessions,
    load_prot_taxid_subset,
    reformat_record,
    stream_reformatted_fasta,
)


def test_reformat_record_emits_bare_taxid_header():
    # Kaiju's mkbwt expects `>TAXID` (a single integer); any other
    # format is mis-parsed as the trailing numeric portion of the
    # accession token. Verified empirically with a synthetic read
    # that round-tripped through `kaiju|<taxid>|<acc>` returning
    # the accession digits rather than the taxid.
    assert reformat_record(">YP_009144834.1 hypothetical protein", 100) == ">100"
    assert reformat_record(">NP_001234.5 description text", 42) == ">42"


def test_collect_wanted_accessions_indexes_versioned_and_base(tmp_path):
    faa = tmp_path / "proteins.faa"
    faa.write_text(
        ">YP_001.1 protein one\n"
        "MAA\n"
        ">YP_002.1 protein two\n"
        "MBB\n"
    )
    wanted = collect_wanted_accessions(faa)
    assert wanted == {"YP_001.1", "YP_001", "YP_002.1", "YP_002"}


def test_load_prot_taxid_subset_filters_to_wanted(tmp_path):
    gz = tmp_path / "prot.accession2taxid.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("accession\taccession.version\ttaxid\tgi\n")
        fh.write("YP_001\tYP_001.1\t100\t111\n")
        fh.write("YP_002\tYP_002.1\t200\t222\n")
        fh.write("YP_999\tYP_999.1\t999\t999\n")  # not wanted

    out = load_prot_taxid_subset(gz, {"YP_001", "YP_001.1", "YP_002"})
    assert out["YP_001"] == 100
    assert out["YP_001.1"] == 100
    assert out["YP_002"] == 200
    assert "YP_999" not in out


def test_stream_reformatted_fasta_emits_kaiju_headers(tmp_path):
    faa = tmp_path / "proteins.faa"
    faa.write_text(
        ">YP_001.1 known protein\n"
        "MAA\n"
        "ABC\n"
        ">YP_002.1 known too\n"
        "MBB\n"
        ">YP_unmapped.1 not in map\n"
        "MZZ\n"
    )
    out = tmp_path / "reformatted.faa"
    written, dropped = stream_reformatted_fasta(
        faa, {"YP_001.1": 100, "YP_001": 100, "YP_002.1": 200}, out
    )
    assert written == 2
    assert dropped == 1
    text = out.read_text()
    assert ">100\n" in text
    assert ">200\n" in text
    # Sequence lines for kept records survive in order.
    assert "MAA" in text
    assert "ABC" in text
    assert "MBB" in text
    # Unmapped record's header and sequence are excluded.
    assert "YP_unmapped" not in text
    assert "MZZ" not in text
