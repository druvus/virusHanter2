"""Tests for the machine-readable provenance sidecar builder.

``run_provenance_<batch>.json`` is the contract the HTML report renders:
databases (short path + robust build identity + date) and the resolved
tool versions that actually ran. It must never carry an absolute path.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.write_provenance import build_provenance  # noqa: E402


def _cfg(tmp_path: Path) -> dict:
    checkv = tmp_path / "checkv" / "checkv-db-v1.5"
    checkv.mkdir(parents=True)
    return {"CHECKV_DB": str(checkv), "HOST_REMOVAL": "bwa"}


def test_build_provenance_shape(tmp_path):
    software_rows = [
        {"env": "fastp", "package": "fastp", "version": "0.24.0", "build": "b0"},
        {"env": "kraken", "package": "kraken2", "version": "2.1.3", "build": "b1"},
        {"env": "kraken", "package": "python", "version": "3.12.2", "build": "h1"},
    ]
    prov = build_provenance(
        _cfg(tmp_path),
        software_rows,
        run_name="20260101_batch",
        assemblers=["MEGAHIT", "metaSPAdes"],
        reporthanter_version="0.9.0",
        snakemake_version="9.23.1",
        python_version="3.12.2",
        generated_utc="2026-07-01T15:00:00Z",
    )

    assert prov["run_name"] == "20260101_batch"
    assert prov["generated_utc"] == "2026-07-01T15:00:00Z"
    assert prov["host_removal_tool"] == "bwa"
    assert prov["assemblers_used"] == ["MEGAHIT", "metaSPAdes"]
    assert prov["reporthanter_version"] == "0.9.0"

    # Databases carry short path + identity.
    checkv = next(d for d in prov["databases"] if d["key"] == "CHECKV_DB")
    assert checkv["path"] == "checkv/checkv-db-v1.5"
    assert checkv["identity"] == "checkv-db-v1.5"

    # Headline map is the tools, not incidental deps.
    assert prov["software_headline"] == {"fastp": "0.24.0", "kraken2": "2.1.3"}
    # Full software table is preserved.
    assert any(r["package"] == "python" for r in prov["software"])


def test_build_provenance_never_leaks_absolute_path(tmp_path):
    prov = build_provenance(
        _cfg(tmp_path),
        [],
        run_name="b",
        assemblers=["MEGAHIT"],
        reporthanter_version="",
        snakemake_version="",
        python_version="",
        generated_utc="t",
    )
    import json

    blob = json.dumps(prov)
    assert str(tmp_path) not in blob
