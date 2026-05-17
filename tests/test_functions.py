"""Regression tests for the sample-discovery helpers in scripts/functions.py.

These cover the gzipped-FASTQ case that the previous, $-anchored regex
silently dropped. Run with `pytest tests/` from the repository root.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.functions import common_suffix, paired_reads  # noqa: E402


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
