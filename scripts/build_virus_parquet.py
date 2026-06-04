"""
Build the virus parquet consumed by `bwa_align_to_kraken_hits`.

The parquet has three columns: `name`, `sequence`, `tax_id`.
`bwa_align_to_kraken_hits` joins the per-classifier viral tax_ids
against this table to pick reference sequences to map reads against.

Two source modes are supported. Default is NCBI Virus, which is
broader than viral RefSeq and includes the ICTV reclassifications
that show up in Kraken / Kaiju but used to be absent from the
parquet. The historical RefSeq mode is preserved for reproducing
older builds.

  --source refseq      Concatenated viral RefSeq FASTA (the
                       pre-2026-05-21 behaviour).
  --source ncbi-virus  All-viral GenBank nucleotide FASTA from
                       NCBI Virus (default). Larger upstream input;
                       paired with --one-rep-per-taxid the output
                       parquet stays in the same size range as the
                       RefSeq build but with much wider taxid
                       coverage.

By default the builder keeps **one longest sequence per tax_id**
(`--one-rep-per-taxid`). The "longest" heuristic biases toward
complete genomes and keeps the parquet to a few rows-per-thousand
of unique taxids regardless of how many GenBank submissions a
species has.

Inputs:
  --fasta   one or more nucleotide FASTA files. Headers must look
            like ``>ACCESSION[.VERSION] free-form description ...``.
  --taxid   gzipped tab-separated NCBI accession2taxid mapping file
            (e.g. ``nucl_gb.accession2taxid.gz``). Columns:
            accession  accession.version  taxid  gi. If omitted,
            `tax_id` is set to 0 for every row — useful only
            for smoke tests.
  --out     output parquet path.

A `build_stats.json` sidecar is written next to the parquet
recording the source, build date, record count, unique tax_id
count and basic length statistics.
"""

import argparse
import gzip
import json
import logging
import sys
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

# When this script is executed directly (`python scripts/build_virus_parquet.py`)
# the script directory is on sys.path but the project root is not. Insert the
# project root so that `from scripts.functions import ...` resolves correctly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.functions import (  # noqa: E402
    find_genus_taxid,
    find_species_taxid,
    parse_nodes_dmp,
)

# Lazy-import pyfastx inside the FASTA-reading helper so consumers
# that only call the taxdump parser / walk-up helpers (e.g. the
# bwa_align_to_kraken_hits script in an env without pyfastx) do not
# need the dependency.


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--fasta",
        nargs="+",
        required=True,
        type=Path,
        help="One or more nucleotide FASTA files.",
    )
    p.add_argument(
        "--taxid",
        type=Path,
        default=None,
        help="NCBI accession2taxid.gz (e.g. nucl_gb.accession2taxid.gz).",
    )
    p.add_argument(
        "--taxdump-nodes",
        type=Path,
        default=None,
        help=(
            "Optional path to an uncompressed NCBI ``nodes.dmp`` from "
            "``taxdump.tar.gz``. When supplied, the parquet is enriched "
            "with `rank` and `genus_taxid` columns so the runtime "
            "coverage rule can apply a rank filter and a genus walk-up. "
            "When omitted, the parquet keeps the three-column "
            "(name, sequence, tax_id) shape and downstream rules behave "
            "as today."
        ),
    )
    p.add_argument("--out", required=True, type=Path, help="Output parquet path.")
    p.add_argument(
        "--source",
        choices=["refseq", "ncbi-virus"],
        default="ncbi-virus",
        help=(
            "Upstream source the FASTA was drawn from. Used only as "
            "provenance in the build_stats.json sidecar; does not "
            "change the parsing logic."
        ),
    )
    p.add_argument(
        "--one-rep-per-taxid",
        dest="one_rep_per_taxid",
        action="store_true",
        default=True,
        help=(
            "Group by tax_id and keep only the longest sequence per "
            "group (default on). Drops rows with tax_id == 0."
        ),
    )
    p.add_argument(
        "--no-one-rep-per-taxid",
        dest="one_rep_per_taxid",
        action="store_false",
        help=(
            "Disable the one-rep-per-taxid filter and write every "
            "FASTA record. Used to reproduce the pre-2026-05-21 "
            "behaviour."
        ),
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return p.parse_args()


def fasta_records(fastas: list[Path]) -> list[tuple[str, str]]:
    """Read every record from every FASTA in `fastas`. Returns a
    list of (header, sequence) tuples in input order. `header` is
    the full header string from the FASTA without the leading '>'.
    """
    import pyfastx  # noqa: PLC0415 - lazy import; see module-level note

    rows: list[tuple[str, str]] = []
    for f in fastas:
        logging.info("Reading FASTA: %s", f)
        fastx = pyfastx.Fastx(str(f))
        for rec in fastx:
            # pyfastx 2.x yields tuples; field 0 is the full header (no '>'),
            # field 1 is the sequence.
            header, seq = rec[0], rec[1]
            rows.append((header, seq))
    logging.info("Total FASTA records read: %d", len(rows))
    return rows


def _accessions_from_header(header: str) -> tuple[str, str]:
    """Return (accession, base_accession) for a FASTA header.

    The accession is the first whitespace-delimited token, with any
    leading "db|" prefix stripped. The base accession is that token
    with its trailing ".VERSION" suffix removed. NCBI's
    accession2taxid file indexes by both, so probing both lets
    headers that drop the version still resolve.
    """
    token = header.split()[0]
    acc = token.split("|")[-1]
    base = acc.split(".")[0]
    return acc, base


def load_taxid_map_subset(path: Path, wanted: set[str]) -> dict[str, int]:
    """Stream the NCBI accession2taxid.gz file and keep only rows
    whose accession or accession.version appears in `wanted`.

    The full `nucl_gb.accession2taxid.gz` has hundreds of millions
    of rows; the prior "load everything into a dict" approach
    allocates >100 GB on a fresh NCBI release. The viral FASTA we
    care about has at most a few hundred thousand records, so the
    kept subset is small and fits in RAM with room to spare.
    """
    if path is None:
        return {}

    logging.info(
        "Streaming taxid map: %s (filtering to %d wanted accessions)",
        path,
        len(wanted),
    )
    out: dict[str, int] = {}
    seen = 0
    with gzip.open(path, "rt") as fh:
        next(fh)  # header line: accession  accession.version  taxid  gi
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc, acc_ver = parts[0], parts[1]
            if acc in wanted or acc_ver in wanted:
                taxid = int(parts[2])
                out[acc] = taxid
                out[acc_ver] = taxid
            seen += 1
            if seen % 50_000_000 == 0:
                logging.info(
                    "  scanned %d lines, kept %d/%d wanted accessions",
                    seen,
                    len({k for k in out if k in wanted}),
                    len(wanted),
                )
    logging.info(
        "Taxid map subset built: %d entries from %d source lines",
        len(out),
        seen,
    )
    return out


def enrich_with_taxdump(
    df: pd.DataFrame, nodes: dict[int, tuple[int, str]]
) -> pd.DataFrame:
    """Return ``df`` with two extra columns: ``rank`` and ``genus_taxid``.

    Both columns are derived from ``nodes`` (parsed from NCBI
    ``nodes.dmp``). Rows whose tax_id is missing from ``nodes`` get
    ``rank = 'unknown'`` and ``genus_taxid = 0``. Original column
    order is preserved with the new columns appended.
    """
    if df.empty:
        return df.assign(rank="unknown", genus_taxid=0)
    rank_col: list[str] = []
    genus_col: list[int] = []
    for tid in df["tax_id"].astype(int).tolist():
        node = nodes.get(tid)
        rank_col.append(node[1] if node else "unknown")
        # If the row itself is already at rank 'genus', use its own
        # tax_id rather than walking further up.
        if node and node[1] == "genus":
            genus_col.append(tid)
        else:
            genus_col.append(find_genus_taxid(tid, nodes))
    return df.assign(rank=rank_col, genus_taxid=genus_col)


def pick_longest_per_taxid(df: pd.DataFrame) -> pd.DataFrame:
    """Keep one row per `tax_id` — the longest sequence.

    `tax_id == 0` rows are dropped (the sentinel for accessions
    that the accession2taxid mapping could not resolve). Ties on
    sequence length are broken by input order: the first record
    encountered for a taxid wins, which makes the build
    deterministic for a fixed input order.
    """
    if df.empty:
        return df.iloc[0:0].copy()
    work = df.loc[df["tax_id"] != 0].copy()
    if work.empty:
        return work
    work["_seqlen"] = work["sequence"].str.len()
    # Stable sort by length descending; idxmax on a sorted-by-length
    # frame gives the first occurrence of the max, i.e. input-order
    # tie-break.
    work = work.sort_values(
        ["tax_id", "_seqlen"],
        ascending=[True, False],
        kind="mergesort",
    )
    picked = work.drop_duplicates(subset=["tax_id"], keep="first")
    return picked.drop(columns=["_seqlen"]).reset_index(drop=True)


def build(
    fastas: list[Path],
    taxid_path: Path | None,
    *,
    one_rep_per_taxid: bool,
    taxdump_nodes_path: Path | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Construct the parquet DataFrame and the build_stats dict."""
    rows = fasta_records(fastas)

    wanted: set[str] = set()
    for header, _ in rows:
        acc, base = _accessions_from_header(header)
        wanted.add(acc)
        wanted.add(base)

    taxid_map = load_taxid_map_subset(taxid_path, wanted) if taxid_path else {}

    names: list[str] = []
    sequences: list[str] = []
    tax_ids: list[int] = []

    missing = 0
    for header, seq in rows:
        acc, base = _accessions_from_header(header)
        tax = taxid_map.get(acc) or taxid_map.get(base)
        if tax is None:
            tax = 0
            missing += 1
        names.append(header)
        sequences.append(seq)
        tax_ids.append(tax)

    if missing:
        logging.warning(
            "No taxid for %d / %d records (taxid set to 0)", missing, len(rows)
        )

    raw = pd.DataFrame({"name": names, "sequence": sequences, "tax_id": tax_ids})

    if one_rep_per_taxid:
        df = pick_longest_per_taxid(raw)
    else:
        df = raw

    rank_distribution: dict[str, int] = {}
    with_genus_taxid_count = 0
    if taxdump_nodes_path is not None:
        logging.info("Parsing taxdump nodes: %s", taxdump_nodes_path)
        nodes = parse_nodes_dmp(taxdump_nodes_path)
        logging.info("Parsed %d taxdump nodes", len(nodes))
        df = enrich_with_taxdump(df, nodes)
        rank_distribution = (
            df["rank"].value_counts().to_dict() if not df.empty else {}
        )
        # Cast counts to int — value_counts returns numpy.int64 which
        # json.dump cannot serialise.
        rank_distribution = {str(k): int(v) for k, v in rank_distribution.items()}
        with_genus_taxid_count = int((df["genus_taxid"] != 0).sum())

    stats = {
        "input_records": int(len(raw)),
        "output_records": int(len(df)),
        "unique_taxids": int(df["tax_id"].nunique()),
        "unresolved_taxid_records": int(missing),
        "mean_sequence_length": float(df["sequence"].str.len().mean()) if len(df) else 0.0,
        "median_sequence_length": float(df["sequence"].str.len().median()) if len(df) else 0.0,
        "one_rep_per_taxid": bool(one_rep_per_taxid),
        "taxdump_enriched": bool(taxdump_nodes_path is not None),
        "rank_distribution": rank_distribution,
        "with_genus_taxid_count": with_genus_taxid_count,
    }
    return df, stats


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    df, stats = build(
        args.fasta,
        args.taxid,
        one_rep_per_taxid=args.one_rep_per_taxid,
        taxdump_nodes_path=args.taxdump_nodes,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logging.info("Wrote %s (%d rows)", args.out, len(df))

    stats_path = args.out.with_name(args.out.stem + "_build_stats.json")
    stats.update(
        {
            "source": args.source,
            "build_date_utc": datetime.now(UTC).isoformat(timespec="seconds"),
            "fasta_inputs": [str(p) for p in args.fasta],
            "taxid_input": str(args.taxid) if args.taxid else None,
            "output_parquet": str(args.out),
        }
    )
    with open(stats_path, "w") as fh:
        json.dump(stats, fh, indent=2)
    logging.info("Wrote build stats %s", stats_path)


if __name__ == "__main__":
    main()
