"""
Build the virus parquet consumed by `bwa_align_to_kraken_hits`.

The parquet has three columns: `name`, `sequence`, `tax_id`. `bwa_align_to_kraken_hits`
joins the Kraken top-20 viral tax_ids against this table to pick reference
sequences to map reads against.

Inputs:
  --fasta   one or more FASTA files (e.g. refreshed NCBI RefSeq viral
            release). Headers must look like
            ``>ACCESSION[.VERSION] free-form description ...``.
  --taxid   gzipped tab-separated NCBI accession2taxid mapping file (e.g.
            `nucl_gb.accession2taxid.gz`). Columns:
            accession  accession.version  taxid  gi
            If omitted, `tax_id` is set to 0 for every row — useful only
            for smoke tests.
  --out     output parquet path.

Re-run this script after every viral RefSeq refresh.
"""

import argparse
import gzip
import logging
from pathlib import Path

import pandas as pd
import pyfastx


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", nargs="+", required=True, type=Path,
                   help="One or more nucleotide FASTA files.")
    p.add_argument("--taxid", type=Path, default=None,
                   help="NCBI accession2taxid.gz (e.g. nucl_gb.accession2taxid.gz).")
    p.add_argument("--out", required=True, type=Path,
                   help="Output parquet path.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args()


def fasta_records(fastas: list[Path]) -> list[tuple[str, str]]:
    """Read every record from every FASTA in `fastas`. Returns a list of
    (header, sequence) tuples in input order. `header` is the full header
    string from the FASTA without the leading '>'.
    """
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
    leading "db|" prefix stripped. The base accession is that token with
    its trailing ".VERSION" suffix removed. NCBI's accession2taxid file
    indexes by both, so probing both lets headers that drop the version
    still resolve.
    """
    token = header.split()[0]
    acc = token.split("|")[-1]
    base = acc.split(".")[0]
    return acc, base


def load_taxid_map_subset(
    path: Path, wanted: set[str]
) -> dict[str, int]:
    """Stream the NCBI accession2taxid.gz file and keep only rows whose
    accession or accession.version appears in `wanted`.

    The full `nucl_gb.accession2taxid.gz` has ~500 million rows; the prior
    "load everything into a dict" approach allocates >100 GB on a fresh
    NCBI release. The viral FASTA we care about has at most ~30k records,
    so the kept subset is small and fits in RAM with room to spare.
    """
    if path is None:
        return {}

    logging.info("Streaming taxid map: %s (filtering to %d wanted accessions)",
                 path, len(wanted))
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
                    seen, len({k for k in out if k in wanted}), len(wanted),
                )
    logging.info("Taxid map subset built: %d entries from %d source lines",
                 len(out), seen)
    return out


def build(fastas: list[Path], taxid_path: Path | None) -> pd.DataFrame:
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
        logging.warning("No taxid for %d / %d records (taxid set to 0)",
                        missing, len(rows))

    return pd.DataFrame({"name": names, "sequence": sequences, "tax_id": tax_ids})


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s [%(levelname)s] %(message)s")

    df = build(args.fasta, args.taxid)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logging.info("Wrote %s (%d rows)", args.out, len(df))


if __name__ == "__main__":
    main()
