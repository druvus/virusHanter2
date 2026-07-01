"""Reference-database provenance helpers (build identity + short paths).

Shared by ``aggregate_run_information`` (run-info CSV columns) and
``write_provenance`` (the machine-readable sidecar the report renders).
Standard library only, so it imports cleanly from any per-rule conda env.

Two concerns:

* **Build identity.** State which DB snapshot produced a result with a
  robust stamp rather than a fragile file mtime. Preference order per
  DB: an explicit ``build_stats.json`` sidecar (written by the refresh
  workflow, carrying ``build_date_utc`` + ``source``), then the
  version-bearing directory / prefix basename (e.g. ``checkv-db-v1.5``,
  ``k2_pluspf_20240112``), with a representative-file mtime as the date
  of last resort.
* **Short paths.** Never surface the operator's absolute filesystem
  layout: only the parent folder and leaf are shown (``checkv/checkv-db-v1.5``).
"""

from __future__ import annotations

import datetime as _dt
import json as _json
from pathlib import Path

# Config keys whose values name a reference database the workflow used.
DB_CONFIG_KEYS = (
    "HUMAN_INDEX",
    "KAIJU_DB",
    "KRAKEN_DB",
    "BLASTN_DB",
    "CHECKV_DB",
    "VIRUS_PARQUET",
    "GENOMAD_DB",
    "TAXDUMP_NODES",
    "SECONDARY_HOST_INDEX",
)

# Directory-backed DBs: representative inner file whose mtime reflects the
# last rebuild (rather than a directory touch).
_DIR_REP_FILES: dict[str, list[str]] = {
    "KAIJU_DB": ["*.fmi", "nodes.dmp"],
    "KRAKEN_DB": ["taxo.k2d", "hash.k2d"],
    "CHECKV_DB": ["genome_db/checkv_reps.dmnd", "checkv-db-v*.tsv"],
    "GENOMAD_DB": ["names.dmp", "genomad_db.dmnd"],
}

# Prefix-backed DBs (BWA / BLAST): sibling suffixes to probe for mtime.
_PREFIX_SUFFIXES: dict[str, list[str]] = {
    "HUMAN_INDEX": [".bwt", ".0123"],
    "BLASTN_DB": [".nhr", ".nal"],
    "SECONDARY_HOST_INDEX": [".bwt", ".0123"],
}


def short_path(path: str) -> str:
    """Return only the parent folder and leaf of ``path``.

    ``/data/refdbs/checkv/checkv-db-v1.5`` -> ``checkv/checkv-db-v1.5``.
    A bare name is returned unchanged; an empty value yields ``""``.
    """
    if not path:
        return ""
    p = Path(path)
    parts = p.parts
    if len(parts) <= 1:
        return p.name or path
    return f"{parts[-2]}/{parts[-1]}"


def _mtime_date(p: Path) -> str:
    try:
        ts = p.stat().st_mtime
        return _dt.datetime.fromtimestamp(ts, _dt.UTC).date().isoformat()
    except OSError:
        return ""


def _read_build_stats(path: Path) -> dict | None:
    try:
        with path.open() as fh:
            data = _json.load(fh)
        return data if isinstance(data, dict) else None
    except (OSError, ValueError):
        return None


def _find_dir_build_stats(root: Path) -> dict | None:
    """Return the first ``*build_stats.json`` in ``root`` (refresh stamp)."""
    for pattern in ("build_stats.json", "*_build_stats.json"):
        for candidate in sorted(root.glob(pattern)):
            data = _read_build_stats(candidate)
            if data:
                return data
    return None


def _stamp_identity(data: dict, fallback: str) -> tuple[str, str]:
    """Render (date, identity) from a build_stats mapping."""
    date = str(data.get("build_date_utc", ""))[:10]
    source = str(data.get("source", "")).strip()
    records = data.get("output_records") or data.get("input_records")
    bits = [b for b in (source, date) if b]
    if records:
        bits.append(f"{records} refs")
    identity = " ".join(bits).strip() or fallback
    return date, identity


def _parquet_identity(prefix: str) -> dict[str, str]:
    parquet = Path(prefix)
    sidecar = parquet.with_name(parquet.stem + "_build_stats.json")
    entry = {"path": short_path(prefix), "date": "", "identity": parquet.name}
    data = _read_build_stats(sidecar) if sidecar.is_file() else None
    if data:
        entry["date"], entry["identity"] = _stamp_identity(data, parquet.name)
    if not entry["date"] and parquet.is_file():
        entry["date"] = _mtime_date(parquet)
    return entry


def _dir_identity(key: str, root: str) -> dict[str, str] | None:
    rp = Path(root)
    if not rp.is_dir():
        return None
    entry = {"path": short_path(root), "date": "", "identity": rp.name}
    data = _find_dir_build_stats(rp)
    if data:
        entry["date"], entry["identity"] = _stamp_identity(data, rp.name)
    if not entry["date"]:
        picked: Path | None = None
        for pattern in _DIR_REP_FILES.get(key, []):
            matches = sorted(rp.glob(pattern))
            if matches:
                picked = matches[0]
                break
        entry["date"] = _mtime_date(picked if picked is not None else rp)
    return entry


def _prefix_identity(key: str, prefix: str) -> dict[str, str]:
    entry = {"path": short_path(prefix), "date": "", "identity": Path(prefix).name}
    for suffix in _PREFIX_SUFFIXES.get(key, []):
        probe = Path(prefix + suffix)
        if probe.exists():
            entry["date"] = _mtime_date(probe)
            break
    return entry


def db_build_identity(cfg: dict) -> dict[str, dict[str, str]]:
    """Resolve a per-DB ``{path, date, identity}`` mapping for each
    configured reference. Unset / unresolvable DBs are dropped.
    """
    out: dict[str, dict[str, str]] = {}

    parquet = cfg.get("VIRUS_PARQUET", "")
    if parquet:
        out["VIRUS_PARQUET"] = _parquet_identity(parquet)

    for key in ("KAIJU_DB", "KRAKEN_DB", "CHECKV_DB", "GENOMAD_DB"):
        root = cfg.get(key, "")
        if not root:
            continue
        entry = _dir_identity(key, root)
        if entry is not None:
            out[key] = entry

    for key in ("HUMAN_INDEX", "BLASTN_DB", "SECONDARY_HOST_INDEX"):
        prefix = cfg.get(key, "")
        if prefix:
            out[key] = _prefix_identity(key, prefix)

    taxdump = cfg.get("TAXDUMP_NODES", "")
    if taxdump:
        path = Path(taxdump)
        if path.is_file():
            out["TAXDUMP_NODES"] = {
                "path": short_path(taxdump),
                "date": _mtime_date(path),
                "identity": short_path(taxdump),
            }

    return out


def databases_used_string(identity: dict[str, dict[str, str]]) -> str:
    """``KEY=folder/leaf`` per DB, semicolon-delimited (short paths only)."""
    return ";".join(f"{k}={v['path']}" for k, v in identity.items() if v.get("path"))


def databases_provenance_string(identity: dict[str, dict[str, str]]) -> str:
    """``KEY=YYYY-MM-DD`` per DB, semicolon-delimited; empty dates dropped."""
    return ";".join(f"{k}={v['date']}" for k, v in identity.items() if v.get("date"))


def databases_build_identity_string(identity: dict[str, dict[str, str]]) -> str:
    """``KEY=<identity>`` per DB, semicolon-delimited; empty ids dropped."""
    return ";".join(f"{k}={v['identity']}" for k, v in identity.items() if v.get("identity"))


def databases_provenance_span_days(identity: dict[str, dict[str, str]]) -> int:
    """Span in days between the oldest and newest resolved DB date.

    Returns 0 when fewer than two dates resolve. A large span means the
    classifier DBs likely came from divergent snapshots.
    """
    parsed: list[_dt.date] = []
    for entry in identity.values():
        value = entry.get("date", "")
        try:
            parsed.append(_dt.date.fromisoformat(value))
        except ValueError:
            continue
    if len(parsed) < 2:
        return 0
    return (max(parsed) - min(parsed)).days
