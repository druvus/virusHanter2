"""Tests for database build-identity resolution and short-path display.

Provenance must state which reference DB snapshot produced a result with
a robust build identity (not a fragile file mtime), and must never leak
the operator's absolute filesystem layout -- only the parent folder and
leaf are shown.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.provenance import (  # noqa: E402
    databases_build_identity_string,
    databases_provenance_span_days,
    databases_used_string,
    db_build_identity,
    short_path,
)


def test_short_path_keeps_only_folder_and_leaf():
    assert short_path("/data/refdbs/checkv/checkv-db-v1.5") == "checkv/checkv-db-v1.5"
    assert short_path("/data/refdbs/human/human_gencode") == "human/human_gencode"


def test_short_path_handles_bare_and_empty():
    assert short_path("") == ""
    assert short_path("all_viruses.parquet") == "all_viruses.parquet"


def test_parquet_identity_prefers_build_stats(tmp_path):
    parquet = tmp_path / "all_viruses.parquet"
    parquet.write_text("x")
    sidecar = tmp_path / "all_viruses_build_stats.json"
    sidecar.write_text(
        json.dumps(
            {"build_date_utc": "2026-05-17T10:00:00", "source": "refseq", "output_records": 14899}
        )
    )
    ident = db_build_identity({"VIRUS_PARQUET": str(parquet)})
    assert ident["VIRUS_PARQUET"]["date"] == "2026-05-17"
    # Identity carries the meaningful build stamp, not just the mtime.
    assert "refseq" in ident["VIRUS_PARQUET"]["identity"]
    assert "14899" in ident["VIRUS_PARQUET"]["identity"]
    # Path is shortened.
    assert ident["VIRUS_PARQUET"]["path"] == "".join(
        [tmp_path.name, "/", "all_viruses.parquet"]
    )


def test_checkv_identity_from_directory_name(tmp_path):
    checkv = tmp_path / "checkv-db-v1.5"
    checkv.mkdir()
    (checkv / "genome_db").mkdir()
    (checkv / "genome_db" / "checkv_reps.dmnd").write_text("x")
    ident = db_build_identity({"CHECKV_DB": str(checkv)})
    # The version-bearing directory name is the identity.
    assert ident["CHECKV_DB"]["identity"] == "checkv-db-v1.5"


def test_kraken_identity_from_directory_name(tmp_path):
    kraken = tmp_path / "k2_pluspf_20240112"
    kraken.mkdir()
    (kraken / "hash.k2d").write_text("x")
    (kraken / "taxo.k2d").write_text("x")
    ident = db_build_identity({"KRAKEN_DB": str(kraken)})
    assert ident["KRAKEN_DB"]["identity"] == "k2_pluspf_20240112"


def test_mtime_fallback_when_no_stamp(tmp_path):
    kaiju = tmp_path / "kaiju_refseq_viral"
    kaiju.mkdir()
    fmi = kaiju / "kaiju_db.fmi"
    fmi.write_text("x")
    old = time.time() - 86400
    os.utime(fmi, (old, old))
    ident = db_build_identity({"KAIJU_DB": str(kaiju)})
    # A date is still resolved from the representative file mtime.
    assert ident["KAIJU_DB"]["date"]


def test_compact_strings_use_short_paths_and_drop_full(tmp_path):
    checkv = tmp_path / "checkv-db-v1.5"
    checkv.mkdir()
    ident = db_build_identity({"CHECKV_DB": str(checkv)})
    used = databases_used_string(ident)
    build = databases_build_identity_string(ident)
    # No absolute path leaks into either cell.
    assert str(tmp_path) not in used
    assert str(tmp_path) not in build
    assert "CHECKV_DB=checkv-db-v1.5" in build


def test_span_days_zero_with_single_db(tmp_path):
    checkv = tmp_path / "checkv-db-v1.5"
    checkv.mkdir()
    ident = db_build_identity({"CHECKV_DB": str(checkv)})
    assert databases_provenance_span_days(ident) == 0
