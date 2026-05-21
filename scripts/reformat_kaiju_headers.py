"""Rewrite a viral protein FASTA's headers for Kaiju's index builder.

Kaiju's ``kaiju-mkbwt`` expects each input record to carry the NCBI
taxonomy id in the header, in the shape ``>kaiju|<taxid>|<accession>``.
The RefSeq viral protein FASTAs as published by NCBI use the regular
``>YP_009144834.1 protein description ...`` shape, so we reformat the
stream before handing it to ``kaiju-mkbwt``.

The taxid for each protein accession comes from NCBI's
``prot.accession2taxid.gz`` (columns: accession, accession.version,
taxid, gi). The file has hundreds of millions of rows so we stream it
once, filtering to only the accessions present in the input FASTA;
this matches the bounded-memory pattern used by
``scripts/build_virus_parquet.load_taxid_map_subset``.

Records whose accession does not resolve to a taxid are dropped from
the output and counted; the unresolved count is logged to stderr.

Usage:

    python scripts/reformat_kaiju_headers.py \\
        --proteins refseq_viral_proteins.faa \\
        --prot-taxid prot.accession2taxid.gz \\
        --out kaiju_refseq_viral.faa
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path


def _accession_from_header(header: str) -> tuple[str, str]:
    """Return (accession, base_accession) for a FASTA header line.

    The accession is the first whitespace-delimited token, with any
    leading ``db|`` prefix stripped. The base accession is that
    token with its trailing ``.VERSION`` suffix removed. NCBI's
    ``prot.accession2taxid`` indexes by both, so probing both lets
    headers that drop the version still resolve.
    """
    token = header.lstrip(">").split()[0]
    acc = token.split("|")[-1]
    base = acc.split(".")[0]
    return acc, base


def reformat_record(header: str, taxid: int) -> str:
    """Return the Kaiju-style replacement header for a record.

    ``header`` is the raw FASTA header line (with or without the
    leading ``>``). ``taxid`` is the NCBI taxonomy id resolved for
    the record's accession.

    Kaiju's ``kaiju-mkbwt`` parses the FASTA header as a single
    integer taxid — anything else (including the
    ``kaiju|<taxid>|<accession>`` format documented on some
    third-party pages) is mis-parsed as the trailing numeric
    portion of the accession. The original accession is dropped;
    Kaiju doesn't need it because classification reports only the
    taxid.
    """
    return f">{taxid}"


def collect_wanted_accessions(faa_path: Path) -> set[str]:
    """Scan the protein FASTA once for the accessions we need to
    resolve. Returns the set of (versioned + unversioned) tokens.

    The scan only reads header lines, so the FASTA stream is cheap
    even at the typical viral RefSeq protein size of a few hundred
    megabytes.
    """
    wanted: set[str] = set()
    with open(faa_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            acc, base = _accession_from_header(line.rstrip("\n"))
            wanted.add(acc)
            wanted.add(base)
    return wanted


def load_prot_taxid_subset(
    path: Path, wanted: set[str]
) -> dict[str, int]:
    """Stream ``prot.accession2taxid.gz`` and keep only rows whose
    accession or accession.version appears in ``wanted``.

    Mirrors ``scripts/build_virus_parquet.load_taxid_map_subset``
    but operates on the protein file. ``prot.accession2taxid`` has
    around a billion rows worldwide; the wanted set bounds memory.
    """
    out: dict[str, int] = {}
    seen = 0
    with gzip.open(path, "rt") as fh:
        next(fh)  # header: accession  accession.version  taxid  gi
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc, acc_ver = parts[0], parts[1]
            if acc in wanted or acc_ver in wanted:
                try:
                    taxid = int(parts[2])
                except ValueError:
                    continue
                out[acc] = taxid
                out[acc_ver] = taxid
            seen += 1
            if seen % 50_000_000 == 0:
                logging.info(
                    "  scanned %d rows; matched %d wanted accessions",
                    seen,
                    len(out),
                )
    return out


def stream_reformatted_fasta(
    faa_path: Path, taxid_map: dict[str, int], out_path: Path
) -> tuple[int, int]:
    """Stream the FASTA and write the reformatted version.

    Each record's header is replaced with the Kaiju-style header.
    Records whose accession does not resolve drop out of the
    output. Returns ``(written, dropped)`` counts.
    """
    written = 0
    dropped = 0
    skip_record = False
    with open(faa_path) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith(">"):
                acc, base = _accession_from_header(line.rstrip("\n"))
                tid = taxid_map.get(acc) or taxid_map.get(base)
                if tid is None:
                    dropped += 1
                    skip_record = True
                    continue
                skip_record = False
                dst.write(reformat_record(line.rstrip("\n"), tid) + "\n")
                written += 1
            else:
                if not skip_record:
                    dst.write(line)
    return written, dropped


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--proteins", required=True, type=Path)
    p.add_argument("--prot-taxid", required=True, type=Path)
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

    logging.info("Scanning protein FASTA for accessions: %s", args.proteins)
    wanted = collect_wanted_accessions(args.proteins)
    logging.info(
        "Found %d unique accession tokens to resolve", len(wanted)
    )

    logging.info(
        "Streaming protein accession2taxid: %s (filtering to wanted accessions)",
        args.prot_taxid,
    )
    taxid_map = load_prot_taxid_subset(args.prot_taxid, wanted)
    logging.info("Resolved %d taxid mappings", len(taxid_map))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    written, dropped = stream_reformatted_fasta(
        args.proteins, taxid_map, args.out
    )
    logging.info(
        "Wrote %d records to %s; dropped %d unresolved records",
        written,
        args.out,
        dropped,
    )


if __name__ == "__main__":
    main()
