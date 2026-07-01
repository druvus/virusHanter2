"""Tests for the resolved-tool-version collector.

Each per-env probe rule dumps the ``conda-meta`` json filename stems of
the exact conda prefix that ran the work (``name-version-build``). The
collector parses those stems -- robustly, since package names can carry
hyphens (``perl-archive-tar``) -- into a flat ``software_versions.tsv``
and a headline ``{package: version}`` map that the sidecar, the run-info
CSV and the report all reuse.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.collect_software_versions import (  # noqa: E402
    ENV_PRIMARY_PACKAGE,
    collect_versions,
    headline_versions,
    parse_conda_meta_stem,
)


def test_parse_simple_stem():
    assert parse_conda_meta_stem("fastp-0.24.0-h5f740d0_0") == (
        "fastp",
        "0.24.0",
        "h5f740d0_0",
    )


def test_parse_hyphenated_package_name():
    # version and build are the last two hyphen-delimited fields; the
    # remainder (which itself contains hyphens) is the package name.
    assert parse_conda_meta_stem("perl-archive-tar-2.40-pl5321hdfd78af_0") == (
        "perl-archive-tar",
        "2.40",
        "pl5321hdfd78af_0",
    )


def test_parse_malformed_stem_returns_whole_as_name():
    # A stem without two hyphens cannot be split; keep it as the name so
    # nothing is silently dropped.
    assert parse_conda_meta_stem("weirdname") == ("weirdname", "", "")


def _write_probe(dirpath: Path, env: str, stems: list[str]) -> Path:
    path = dirpath / f"{env}.tsv"
    path.write_text("".join(f"{env}\t{stem}\n" for stem in stems))
    return path


def test_collect_versions_merges_and_parses(tmp_path):
    p1 = _write_probe(tmp_path, "kraken", ["kraken2-2.1.3-pl5321hdcf5f25_0", "python-3.12.2-h1"])
    p2 = _write_probe(tmp_path, "fastp", ["fastp-0.24.0-h5f740d0_0"])

    rows = collect_versions([p1, p2])

    # Flat rows carry env, package, version, build.
    kraken2 = next(r for r in rows if r["package"] == "kraken2")
    assert kraken2 == {
        "env": "kraken",
        "package": "kraken2",
        "version": "2.1.3",
        "build": "pl5321hdcf5f25_0",
    }
    assert any(r["package"] == "fastp" and r["version"] == "0.24.0" for r in rows)


def test_headline_picks_primary_package_per_env(tmp_path):
    p1 = _write_probe(tmp_path, "kraken", ["kraken2-2.1.3-b0", "python-3.12.2-h1"])
    p2 = _write_probe(tmp_path, "blastn", ["blast-2.16.0-h1", "python-3.12.2-h1"])
    rows = collect_versions([p1, p2])

    headline = headline_versions(rows)

    # The kraken env's headline tool is kraken2 (not python); blastn -> blast.
    assert headline["kraken2"] == "2.1.3"
    assert headline["blast"] == "2.16.0"
    # Incidental packages (python) are not headline tools.
    assert "python" not in headline


def test_primary_package_map_covers_probed_envs():
    # Guard against a probe env being added without a headline mapping.
    for env in ("fastp", "bwa", "samtools", "kraken", "kaiju", "megahit",
                "spades", "pilon", "blastn", "checkv", "mosdepth",
                "multiqc", "genomad", "quast", "hostile"):
        assert env in ENV_PRIMARY_PACKAGE
