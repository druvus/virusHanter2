"""Regression tests for the sample-discovery and Kraken helpers in
``scripts/functions.py``.

These cover the gzipped-FASTQ case that the previous, $-anchored regex
silently dropped, and the D / R1 tax_lvl-anchor cases for the two
Kraken2 DB shapes the pipeline supports (pluspf vs. viral-only). Run
with `pytest tests/` from the repository root.
"""
import io
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.functions import (  # noqa: E402
    common_suffix,
    paired_reads,
    wrangle_kraken,
)


def _touch_all(directory: Path, names: list[str]) -> None:
    for name in names:
        (directory / name).touch()


def test_paired_reads_finds_gzipped(tmp_path):
    _touch_all(tmp_path, [
        "a_R1_001.fastq.gz",
        "a_R2_001.fastq.gz",
        "b_R1_001.fastq.gz",
        "b_R2_001.fastq.gz",
    ])
    assert paired_reads(str(tmp_path)) == ["a_R", "b_R"]


def test_paired_reads_finds_uncompressed(tmp_path):
    _touch_all(tmp_path, [
        "x_R1.fastq",
        "x_R2.fastq",
    ])
    assert paired_reads(str(tmp_path)) == ["x_R"]


def test_common_suffix_handles_gz(tmp_path):
    _touch_all(tmp_path, [
        "a_R1_001.fastq.gz",
        "a_R2_001.fastq.gz",
    ])
    # Differs only in the R{1,2} character; everything from `_001.fastq.gz`
    # to end of name is shared.
    assert common_suffix(str(tmp_path)) == "_001.fastq.gz"


def test_common_suffix_uncompressed(tmp_path):
    _touch_all(tmp_path, [
        "x_R1.fastq",
        "x_R2.fastq",
    ])
    assert common_suffix(str(tmp_path)) == ".fastq"


def test_helpers_skip_non_fastq_files(tmp_path):
    # Mix in a stray .txt and a hidden .DS_Store; they must not break
    # discovery or contaminate the sample list.
    _touch_all(tmp_path, [
        "a_R1.fastq.gz",
        "a_R2.fastq.gz",
        ".DS_Store",
        "notes.txt",
    ])
    assert paired_reads(str(tmp_path)) == ["a_R"]
    assert common_suffix(str(tmp_path)) == ".fastq.gz"


# ---------------------------------------------------------------------------
# wrangle_kraken: D / R1 anchor coverage
#
# Kraken2's pluspf DB labels the "Viruses" row with tax_lvl == "D"
# (Domain). The smaller viral-only DBs (e.g. k2_viral_*) drop it to
# tax_lvl == "R1" because there is no sibling superkingdom. The
# `domain` column carried down by wrangle_kraken must read "Viruses"
# in both cases so the downstream `domain == "Viruses"` filter in
# bwa_align_to_kraken_hits and per_virus_metrics works for either DB.
# ---------------------------------------------------------------------------


def _write_kraken_tsv(path: Path, rows: list[tuple[float, int, int, str, int, str]]) -> None:
    """Write a Kraken2-style report TSV (no header) to ``path``.

    Columns: percent, count_clades, count, tax_lvl, taxonomy_id, name.
    """
    lines = ["\t".join(str(c) for c in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def test_wrangle_kraken_pluspf_shape_domain_is_viruses(tmp_path):
    # pluspf-style: root (R) -> superkingdoms at tax_lvl 'D'.
    tsv = tmp_path / "pluspf.report"
    _write_kraken_tsv(tsv, [
        (5.0, 50, 50, "U", 0, "unclassified"),
        (95.0, 950, 0, "R", 1, "root"),
        (10.0, 100, 0, "D", 2, "Bacteria"),
        (85.0, 850, 0, "D", 10239, "Viruses"),
        (85.0, 850, 50, "S", 10335, "Human alphaherpesvirus 3"),
    ])
    df = wrangle_kraken(str(tsv))
    by_name = dict(zip(df["name"], df["domain"]))
    assert by_name["Viruses"] == "Viruses"
    assert by_name["Human alphaherpesvirus 3"] == "Viruses"
    assert by_name["Bacteria"] == "Bacteria"


def test_wrangle_kraken_viral_only_shape_domain_is_viruses(tmp_path):
    # k2_viral-style: no sibling superkingdom, so Viruses sits at
    # tax_lvl == 'R1' rather than 'D'. The wrangler must still
    # propagate "Viruses" as the domain to every species row below.
    tsv = tmp_path / "k2_viral.report"
    _write_kraken_tsv(tsv, [
        (3.0, 30, 30, "U", 0, "unclassified"),
        (97.0, 970, 0, "R", 1, "root"),
        (97.0, 970, 0, "R1", 10239, "Viruses"),
        (97.0, 970, 0, "R2", 2731341, "Duplodnaviria"),
        (97.0, 970, 50, "S", 10335, "Human alphaherpesvirus 3"),
    ])
    df = wrangle_kraken(str(tsv))
    by_name = dict(zip(df["name"], df["domain"]))
    assert by_name["Viruses"] == "Viruses"
    assert by_name["Duplodnaviria"] == "Viruses"
    assert by_name["Human alphaherpesvirus 3"] == "Viruses"


def test_wrangle_kraken_pluspf_cellular_r1_does_not_shadow_d(tmp_path):
    # In pluspf the only R1 row is "cellular organisms". It briefly
    # carries that name as `domain` until the next D row (Bacteria)
    # overrides; the parity invariant is that every D-and-below row
    # keeps the correct superkingdom as its domain.
    tsv = tmp_path / "pluspf_cellular.report"
    _write_kraken_tsv(tsv, [
        (3.0, 30, 30, "U", 0, "unclassified"),
        (97.0, 970, 0, "R", 1, "root"),
        (50.0, 500, 0, "R1", 131567, "cellular organisms"),
        (40.0, 400, 0, "D", 2, "Bacteria"),
        (40.0, 400, 50, "S", 562, "Escherichia coli"),
        (10.0, 100, 0, "D", 10239, "Viruses"),
        (10.0, 100, 50, "S", 10335, "Human alphaherpesvirus 3"),
    ])
    df = wrangle_kraken(str(tsv))
    by_name = dict(zip(df["name"], df["domain"]))
    assert by_name["Bacteria"] == "Bacteria"
    assert by_name["Escherichia coli"] == "Bacteria"
    assert by_name["Viruses"] == "Viruses"
    assert by_name["Human alphaherpesvirus 3"] == "Viruses"
