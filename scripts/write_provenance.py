"""Write the per-run provenance sidecar (``run_provenance_<batch>.json``).

The sidecar is the single machine-readable record of what produced a
run: the reference databases (short path + robust build identity + date)
and the conda-resolved versions of the tools that actually ran. The HTML
report renders it (via reporthanter's ``--provenance_file``) and it is
diffable across runs. Standard library only.

Per project convention no absolute path is ever written -- database
paths are reduced to their parent folder and leaf by
``scripts.provenance.short_path``.

Note: no ``from __future__ import annotations`` here. Snakemake prepends
its own boilerplate when materialising a ``script:`` rule, so a future
import would no longer sit at the top of the file and would raise a
SyntaxError. The py>=3.11 rule env evaluates the annotations natively.
"""

import csv
import datetime as _dt
import json
from pathlib import Path

from scripts.collect_software_versions import headline_versions, read_versions_table
from scripts.provenance import db_build_identity

# Order databases in the sidecar the way an operator reads a run: host,
# then the classifiers, then the assembly/coverage references.
_DB_ORDER = (
    "HUMAN_INDEX",
    "SECONDARY_HOST_INDEX",
    "KRAKEN_DB",
    "KAIJU_DB",
    "BLASTN_DB",
    "CHECKV_DB",
    "GENOMAD_DB",
    "VIRUS_PARQUET",
    "TAXDUMP_NODES",
)


def build_provenance(
    cfg: dict,
    software_rows: list[dict[str, str]],
    *,
    run_name: str,
    assemblers: list[str],
    reporthanter_version: str,
    snakemake_version: str,
    python_version: str,
    generated_utc: str,
) -> dict:
    """Assemble the provenance mapping from config + resolved versions."""
    identity = db_build_identity(cfg)
    databases = []
    for key in _DB_ORDER:
        entry = identity.get(key)
        if entry:
            databases.append({"key": key, **entry})

    return {
        "run_name": run_name,
        "generated_utc": generated_utc,
        "host_removal_tool": cfg.get("HOST_REMOVAL", "bwa"),
        "assemblers_used": list(assemblers),
        "reporthanter_version": reporthanter_version,
        "snakemake_version": snakemake_version,
        "python_version": python_version,
        "databases": databases,
        "software": list(software_rows),
        "software_headline": headline_versions(software_rows),
    }


def _write_tsv_companion(prov: dict, out_path: Path) -> None:
    """Flat two-table TSV companion for quick eyeballing / grep."""
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["kind", "name", "value", "detail"])
        for db in prov["databases"]:
            writer.writerow(["database", db["key"], db.get("identity", ""), db.get("date", "")])
        for tool, version in sorted(prov["software_headline"].items()):
            writer.writerow(["software", tool, version, ""])


def main() -> None:
    # snakemake is injected into globals by the Snakemake `script:` runner.
    sm = globals()["snakemake"]  # noqa: F821 (provided at runtime)
    cfg = dict(sm.config)

    software_tsv = getattr(sm.input, "software_versions", None)
    # The input is the merged software_versions.tsv (header + rows), so
    # read it back with read_versions_table -- NOT collect_versions,
    # which parses the raw env<TAB>stem probe dumps.
    software_rows = read_versions_table(Path(software_tsv)) if software_tsv else []

    generated_utc = sm.params.get("generated_utc", "") or (
        _dt.datetime.now(_dt.UTC).replace(microsecond=0).isoformat()
    )
    prov = build_provenance(
        cfg,
        software_rows,
        run_name=str(sm.params.run_name),
        assemblers=list(sm.params.get("assemblers", ["MEGAHIT"])),
        reporthanter_version=str(sm.params.get("reporthanter_version", "")),
        snakemake_version=str(sm.params.get("snakemake_version", "")),
        python_version=str(sm.params.get("python_version", "")),
        generated_utc=generated_utc,
    )

    out_json = Path(sm.output.json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w") as fh:
        json.dump(prov, fh, indent=2)
        fh.write("\n")

    tsv_out = getattr(sm.output, "tsv", None)
    if tsv_out:
        _write_tsv_companion(prov, Path(tsv_out))


if __name__ == "__main__":
    main()
