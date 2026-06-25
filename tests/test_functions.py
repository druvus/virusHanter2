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

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd  # noqa: E402

from scripts.functions import (  # noqa: E402
    common_suffix,
    dummy_contig_sentinel,
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


def test_paired_reads_raises_on_odd_file_count(tmp_path):
    # A missing mate (odd number of FASTQ files) must fail with a clear
    # message naming the problem, not an opaque IndexError from samples[i+1].
    _touch_all(tmp_path, [
        "a_R1.fastq.gz",
        "a_R2.fastq.gz",
        "b_R1.fastq.gz",  # b_R2 missing
    ])
    with pytest.raises(ValueError, match="odd number"):
        paired_reads(str(tmp_path))


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


# --- dummy_contig_sentinel: surface a silent total assembly failure ---
# When an assembler produces no usable contigs it writes a DUMMY_CONTIG
# sentinel. CheckV still reports it but BLASTN drops it (no DB hit), so
# the BLAST/CheckV inner join is empty and per_virus_metrics cannot tell
# a silent failure apart from a real negative. dummy_contig_sentinel
# carries one dummy-named row through so the downstream flag can fire.

_MERGED_COLS = ["name", "assembler", "match_name", "tax_id", "viral_genes"]


def test_dummy_contig_sentinel_emits_row_on_failed_assembly():
    merged = pd.DataFrame(columns=_MERGED_COLS)  # empty inner join
    checkv = pd.DataFrame({"name": ["DUMMY_CONTIG_pilon"], "viral_genes": [0]})

    out = dummy_contig_sentinel(merged, checkv, "MEGAHIT")

    assert len(out) == 1
    assert out["name"].iloc[0] == "DUMMY_CONTIG_pilon"
    assert out["assembler"].iloc[0] == "MEGAHIT"
    # Schema is preserved so concatenation downstream stays clean.
    assert list(out.columns) == _MERGED_COLS


def test_dummy_contig_sentinel_noop_when_merge_has_real_hits():
    merged = pd.DataFrame(
        {
            "name": ["NODE_1"],
            "assembler": ["MEGAHIT"],
            "match_name": ["Some virus"],
            "tax_id": [12345],
            "viral_genes": [3],
        }
    )
    checkv = pd.DataFrame({"name": ["NODE_1"], "viral_genes": [3]})

    out = dummy_contig_sentinel(merged, checkv, "MEGAHIT")

    pd.testing.assert_frame_equal(out, merged)


def test_dummy_contig_sentinel_noop_when_checkv_has_real_contigs():
    # A genuine negative: real contigs assembled, none viral, so the
    # join is empty but CheckV lists real (non-dummy) contig names. This
    # must NOT be flagged as an assembly failure.
    merged = pd.DataFrame(columns=_MERGED_COLS)
    checkv = pd.DataFrame({"name": ["NODE_1", "NODE_2"], "viral_genes": [0, 0]})

    out = dummy_contig_sentinel(merged, checkv, "metaSPAdes")

    assert out.empty


def test_dummy_contig_sentinel_noop_when_checkv_empty():
    merged = pd.DataFrame(columns=_MERGED_COLS)
    checkv = pd.DataFrame(columns=["name", "viral_genes"])

    out = dummy_contig_sentinel(merged, checkv, "MEGAHIT")

    assert out.empty
