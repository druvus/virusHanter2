"""Compare VIRUS_PARQUET's taxid universe with a Kraken2 viral DB.

Run ``kraken2-inspect`` against the configured DB to enumerate the
DB's taxa, intersect with the parquet's tax_id column, and emit a
sidecar TSV plus a small JSON summary attached to the existing
``build_stats.json``.

The diagnostic exists because the production Kraken2 viral DB
(``k2_viral_<DATE>``) is built from RefSeq viral; running this on a
RefSeq-built parquet quantifies the residual asymmetry between the
two RefSeq snapshots NCBI publishes for nucleotide vs Kraken2-index
purposes.

Output TSV columns: tax_id, name, in_parquet, in_kraken2, only_in.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_kraken_inspect(text: str) -> dict[int, str]:
    """Parse the standard ``kraken2-inspect`` output into
    ``{taxonomy_id: name}``.

    ``kraken2-inspect`` writes a tab-separated table with six
    columns (matching ``kraken2`` itself): percent, count,
    count_clades, tax_lvl, taxonomy_id, name. The name field is
    indented with two leading spaces per rank-depth; the indentation
    is stripped here.

    Lines starting with ``#`` or that have fewer than six fields are
    ignored.
    """
    out: dict[int, str] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        try:
            tid = int(parts[4].strip())
        except (ValueError, IndexError):
            continue
        name = parts[5].strip()
        out[tid] = name
    return out


def run_kraken_inspect(kraken_db: Path) -> str:
    """Invoke ``kraken2-inspect`` and return its stdout.

    Note: ``--skip-counts`` makes ``kraken2-inspect`` print only the
    summary header (``# Total taxonomy nodes: N``) without the
    per-taxon table, which is the opposite of what we want; the
    full per-taxon listing is the input to :func:`parse_kraken_inspect`.
    """
    result = subprocess.run(
        ["kraken2-inspect", "--db", str(kraken_db)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def build_overlap_frame(
    parquet_tax_ids: set[int],
    kraken_taxid_to_name: dict[int, str],
) -> pd.DataFrame:
    """Build the long-form overlap TSV.

    Includes every taxid in either set, with three boolean
    indicators and a short tag in ``only_in`` for convenience
    (``parquet``, ``kraken2``, or ``both``).
    """
    all_tids = sorted(parquet_tax_ids | set(kraken_taxid_to_name.keys()))
    rows = []
    for tid in all_tids:
        in_p = tid in parquet_tax_ids
        in_k = tid in kraken_taxid_to_name
        only = "both" if in_p and in_k else ("parquet" if in_p else "kraken2")
        rows.append(
            {
                "tax_id": tid,
                "name": kraken_taxid_to_name.get(tid, ""),
                "in_parquet": in_p,
                "in_kraken2": in_k,
                "only_in": only,
            }
        )
    return pd.DataFrame(rows)


def summarise(
    parquet_tax_ids: set[int],
    kraken_taxid_to_name: dict[int, str],
    kraken_db: Path,
) -> dict[str, object]:
    """Return the counters that will be merged into build_stats."""
    kraken_set = set(kraken_taxid_to_name.keys())
    inter = parquet_tax_ids & kraken_set
    parquet_only = parquet_tax_ids - kraken_set
    kraken_only = kraken_set - parquet_tax_ids
    return {
        "kraken2_db_path": str(kraken_db),
        "kraken2_db_taxids": int(len(kraken_set)),
        "parquet_taxids": int(len(parquet_tax_ids)),
        "intersection_count": int(len(inter)),
        "parquet_only_count": int(len(parquet_only)),
        "kraken2_only_count": int(len(kraken_only)),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--parquet", required=True, type=Path)
    p.add_argument("--kraken-db", required=True, type=Path)
    p.add_argument("--build-stats", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stderr,
    )

    logging.info("Loading parquet taxid column: %s", args.parquet)
    df = pd.read_parquet(args.parquet, columns=["tax_id"])
    parquet_tax_ids = set(df["tax_id"].dropna().astype(int).tolist())
    logging.info("Parquet has %d unique taxids", len(parquet_tax_ids))

    logging.info("Running kraken2-inspect on %s", args.kraken_db)
    inspect_text = run_kraken_inspect(args.kraken_db)
    kraken_taxid_to_name = parse_kraken_inspect(inspect_text)
    logging.info("Kraken2 DB has %d taxids", len(kraken_taxid_to_name))

    overlap = build_overlap_frame(parquet_tax_ids, kraken_taxid_to_name)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    overlap.to_csv(args.out, sep="\t", index=False)
    logging.info("Wrote overlap TSV (%d rows) to %s", len(overlap), args.out)

    summary = summarise(parquet_tax_ids, kraken_taxid_to_name, args.kraken_db)
    logging.info(
        "intersection=%d parquet_only=%d kraken2_only=%d",
        summary["intersection_count"],
        summary["parquet_only_count"],
        summary["kraken2_only_count"],
    )

    # Merge into build_stats.json in place.
    if args.build_stats.exists() and args.build_stats.stat().st_size > 0:
        existing = json.loads(args.build_stats.read_text())
    else:
        existing = {}
    existing.update(summary)
    args.build_stats.write_text(json.dumps(existing, indent=2))


if __name__ == "__main__":
    main()
