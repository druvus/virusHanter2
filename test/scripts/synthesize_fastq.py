"""
Synthesize a tiny paired-end FASTQ pair for the virusHanter2 smoke test.

Writes 100 read pairs:
  - 30 derived from the synthetic "host" reference (will map to HUMAN_INDEX)
  - 30 derived from the synthetic viral reference (drive Kraken/Kaiju/BLAST hits)
  - 40 random reads (force unclassified bins to be non-empty)

The output FASTQ files are gzipped so they remain "*.fastq.gz" on disk for
the README convention, even though paired_reads() in scripts/functions.py
ignores the .gz suffix in its regex; the smoke runner generates a sibling
.fastq pair that the workflow can actually pick up.
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

QUAL = "I"  # Phred 40, ASCII 73

HOST_REF = (
    "ACGTGACTGACGTAGCTAGCATCGACTGACTACGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGATCGTAGCT" * 60
)
VIRUS_REF = (
    "TTGCAACGGGCAAATAGTCAGCGGCATTACCTGCAAACGAACAGTATTACCGCAGGCCAGTCTGTAGGAAACAGGGCAAGTGCCATC" * 60
)


def make_read(template: str, start: int, length: int = 150) -> str:
    return template[start : start + length].ljust(length, "N")


def write_pair(
    out_r1: Path,
    out_r2: Path,
    seed: int = 17,
) -> None:
    random.seed(seed)

    records: list[tuple[str, str, str]] = []  # (id, r1, r2)
    # Host reads
    for i in range(30):
        pos = random.randint(0, len(HOST_REF) - 400)
        r1 = make_read(HOST_REF, pos)
        r2 = make_read(HOST_REF, pos + 200)
        records.append((f"host_{i:03d}", r1, r2))
    # Virus reads
    for i in range(30):
        pos = random.randint(0, len(VIRUS_REF) - 400)
        r1 = make_read(VIRUS_REF, pos)
        r2 = make_read(VIRUS_REF, pos + 200)
        records.append((f"virus_{i:03d}", r1, r2))
    # Random reads
    bases = "ACGT"
    for i in range(40):
        r1 = "".join(random.choice(bases) for _ in range(150))
        r2 = "".join(random.choice(bases) for _ in range(150))
        records.append((f"random_{i:03d}", r1, r2))

    random.shuffle(records)

    opener = gzip.open if out_r1.suffix == ".gz" else open
    mode = "wt"
    with opener(out_r1, mode) as f1, opener(out_r2, mode) as f2:
        for read_id, r1, r2 in records:
            f1.write(f"@{read_id}/1\n{r1}\n+\n{QUAL * len(r1)}\n")
            f2.write(f"@{read_id}/2\n{r2}\n+\n{QUAL * len(r2)}\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out-r1", type=Path, required=True)
    p.add_argument("--out-r2", type=Path, required=True)
    args = p.parse_args()

    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r2.parent.mkdir(parents=True, exist_ok=True)
    write_pair(args.out_r1, args.out_r2)


def get_host_reference() -> str:
    return HOST_REF


def get_virus_reference() -> str:
    return VIRUS_REF


if __name__ == "__main__":
    main()
