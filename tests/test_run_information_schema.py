"""Parity guard for the run_information_<batch>.csv column schema.

The first fourteen columns of the run-information CSV are parity-locked
against the original virusHanter (docs/PARITY_NOTES.md, lines 64-67 and
105): their names, order and values must stay byte-identical. A refactor
that reorders or drops one is a silent correctness regression that only a
manual diff would otherwise catch. This test fails loudly instead.

`aggregate_run_information` imports `reporthanter` at module level, which
is not present in the bare base env used for the rest of the unit suite,
so the test skips when reporthanter is unavailable and runs in CI (which
installs it).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("reporthanter")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.aggregate_run_information import aggregate_sample_info  # noqa: E402

# The parity-locked leading schema, in order. Keep identical to
# docs/PARITY_NOTES.md.
PARITY_LOCKED_COLUMNS = [
    "run_name",
    "sample_name",
    "date",
    "read_len",
    "number_reads",
    "mapped_to_human_percent",
    "kraken_virus_percent",
    "kaiju_virus_percent",
    "number_of_contigs",
    "top_contigs_blastn",
    "top_virus_kaiju",
    "html_report",
    "kaiju_report",
    "blastn_report",
]

_FLAGSTAT = """2000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
1500 + 0 mapped (75.00% : N/A)
2000 + 0 paired in sequencing
1000 + 0 read1
1000 + 0 read2
1400 + 0 properly paired (70.00% : N/A)
1400 + 0 with itself and mate mapped
100 + 0 singletons (5.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
"""


def _build_sample(tmp_path: Path) -> Path:
    """Create a minimal but valid sample folder for aggregate_sample_info."""
    run = tmp_path / "20260101_batch"
    sample = run / "s1"
    (sample / "REPORT").mkdir(parents=True)
    (sample / "FASTP").mkdir()
    (sample / "logs").mkdir()
    (sample / "KRAKEN").mkdir()
    (sample / "KAIJU").mkdir()

    (sample / "REPORT" / "s1.html").write_text("<html>report</html>")
    (sample / "FASTP" / "s1.fastp.json").write_text(
        json.dumps(
            {"summary": {"before_filtering": {"read1_mean_length": 150, "total_reads": 2000}}}
        )
    )
    (sample / "logs" / "human_contamination_flagstat.txt").write_text(_FLAGSTAT)
    (sample / "KRAKEN" / "s1.kraken.csv").write_text(
        "percent,count_clades,count,tax_lvl,taxonomy_id,name,domain\n"
        "5.0,50,0,D,10239,Viruses,Viruses\n"
        "5.0,50,50,S,10335,Human alphaherpesvirus 3,Viruses\n"
    )
    (sample / "KAIJU" / "s1.kaiju.table.tsv").write_text(
        "file\tpercent\treads\ttaxon_id\ttaxon_name\n"
        "s1\t2.5\t50\t10335\tHuman alphaherpesvirus 3\n"
    )
    return sample


def test_run_information_leading_columns_are_parity_locked(tmp_path):
    sample = _build_sample(tmp_path)
    df = aggregate_sample_info(sample, assemblers=["MEGAHIT"])

    assert list(df.columns)[: len(PARITY_LOCKED_COLUMNS)] == PARITY_LOCKED_COLUMNS
    # Spot-check parity-locked values resolve from the fixture.
    assert df["run_name"].iloc[0] == "20260101_batch"
    assert df["sample_name"].iloc[0] == "s1"
    assert df["date"].iloc[0] == "20260101"
    assert df["number_reads"].iloc[0] == 2000


def test_run_information_trailing_columns_follow_locked_prefix(tmp_path):
    sample = _build_sample(tmp_path)
    df = aggregate_sample_info(sample, assemblers=["MEGAHIT"])

    # Additive trailing columns must come strictly after the locked 14.
    cols = list(df.columns)
    assert cols[: len(PARITY_LOCKED_COLUMNS)] == PARITY_LOCKED_COLUMNS
    assert "duplicate_pairs" in cols[len(PARITY_LOCKED_COLUMNS) :]
    assert "host_removal_tool" in cols[len(PARITY_LOCKED_COLUMNS) :]
    # Provenance / version columns are additive and follow the locked prefix.
    assert "databases_build_identity" in cols[len(PARITY_LOCKED_COLUMNS) :]
    assert "tool_versions" in cols[len(PARITY_LOCKED_COLUMNS) :]


def test_tool_versions_column_reads_software_versions_tsv(tmp_path):
    sample = _build_sample(tmp_path)
    versions = tmp_path / "software_versions.tsv"
    versions.write_text(
        "env\tpackage\tversion\tbuild\n"
        "fastp\tfastp\t0.24.0\th5f740d0_0\n"
        "kraken\tkraken2\t2.1.3\tpl5321_0\n"
        "kraken\tpython\t3.12.2\th1\n"
    )
    from scripts.aggregate_run_information import _tool_versions_string

    cell = _tool_versions_string(versions)
    # Headline tools only, compact and sorted; incidental python dropped.
    assert cell == "fastp=0.24.0;kraken2=2.1.3"


def test_tool_versions_blank_when_missing(tmp_path):
    from scripts.aggregate_run_information import _tool_versions_string

    assert _tool_versions_string(None) == ""
    assert _tool_versions_string(tmp_path / "absent.tsv") == ""


def test_genomad_columns_populate_from_improved_contigs_path(tmp_path):
    # geNomad names its outputs after the input FASTA stem
    # (<sample>_improved_contigs.fasta), so the rule writes
    # <asm>/GENOMAD/<sample>_improved_contigs_summary/<sample>_improved_contigs_virus_summary.tsv.
    # aggregate_sample_info must stat that exact path; the earlier
    # <sample>_summary/<sample>_virus_summary.tsv form never resolved,
    # leaving the geNomad columns silently blank on GENOMAD: "TRUE" runs.
    sample = _build_sample(tmp_path)
    sample_name = sample.name
    gdir = sample / "MEGAHIT" / "GENOMAD" / f"{sample_name}_improved_contigs_summary"
    gdir.mkdir(parents=True)
    (gdir / f"{sample_name}_improved_contigs_virus_summary.tsv").write_text(
        "seq_name\tvirus_score\nk57_1\t0.95\nk57_2\t0.80\n"
    )

    df = aggregate_sample_info(sample, assemblers=["MEGAHIT"])

    assert df["genomad_viral_contigs"].iloc[0] == 2
    assert df["genomad_max_virus_score"].iloc[0] == 0.95
