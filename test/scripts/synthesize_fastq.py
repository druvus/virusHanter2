"""
Synthesize a tiny paired-end FASTQ pair for the virusHanter2 smoke test.

Writes 100 read pairs:
  - 30 derived from the synthetic "host" reference (will map to HUMAN_INDEX)
  - 30 derived from the synthetic viral reference (drive Kraken/Kaiju/BLAST hits)
  - 40 random reads (force unclassified bins to be non-empty)

The output FASTQ files are gzipped when the requested path ends in
".gz", matching the way production MiSeq batches are stored on disk.
paired_reads() / common_suffix() in scripts/functions.py recognise the
.gz suffix.
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

QUAL = "I"  # Phred 40, ASCII 73


def _pseudo_random_seq(seed: int, length: int) -> str:
    """Deterministic pseudo-random DNA of `length` from `seed`. Used so the
    smoke fixtures are reproducible without committing 10kb of FASTA.
    """
    rng = random.Random(seed)
    bases = "ACGT"
    return "".join(rng.choice(bases) for _ in range(length))


# 5 kb non-repetitive synthetic references. Tandem repeats collapse in the
# de Bruijn graph and MEGAHIT assembles 0 contigs, so the reference must
# carry enough sequence complexity for the assembler to find a path.
HOST_REF = _pseudo_random_seq(seed=1, length=5000)
VIRUS_REF = _pseudo_random_seq(seed=2, length=5000)


def make_read(template: str, start: int, length: int = 150) -> str:
    return template[start : start + length].ljust(length, "N")


def write_pair(
    out_r1: Path,
    out_r2: Path,
    seed: int = 17,
) -> None:
    random.seed(seed)

    # Aim for ~50x coverage of the viral reference so MEGAHIT produces a
    # real contig (avoiding the dummy-contig fallback that confuses CheckV).
    # Host reads stay modest; they only need to drive non-zero mapped counts
    # in the flagstat output.
    records: list[tuple[str, str, str]] = []  # (id, r1, r2)
    # Host reads (60 pairs)
    for i in range(60):
        pos = random.randint(0, len(HOST_REF) - 400)
        r1 = make_read(HOST_REF, pos)
        r2 = make_read(HOST_REF, pos + 200)
        records.append((f"host_{i:03d}", r1, r2))
    # Virus reads (800 pairs ~ 50x coverage at 150bp paired, 5 kb reference)
    for i in range(800):
        pos = random.randint(0, len(VIRUS_REF) - 400)
        r1 = make_read(VIRUS_REF, pos)
        r2 = make_read(VIRUS_REF, pos + 200)
        records.append((f"virus_{i:03d}", r1, r2))
    # Random reads (50 pairs) to keep the unclassified bin non-empty
    bases = "ACGT"
    for i in range(50):
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
