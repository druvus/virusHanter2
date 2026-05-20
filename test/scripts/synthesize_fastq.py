"""
Synthesize a paired-end FASTQ pair for the virusHanter2 smoke test.

A single invocation writes one paired-end sample tagged with a single
viral target. The smoke fixture workflow calls this script once per
sample so the smoke exercises multi-sample paths (per-sample reports,
run-level aggregation, per-virus aggregation).

Three 5 kb non-repetitive synthetic viruses are defined here:

  alpha  taxid 100001
  beta   taxid 100002
  gamma  taxid 100003

Each sample is given one of them; the read mix is
  - small host spike
  - reads from the chosen virus at roughly 50x coverage
  - random reads to keep the unclassified bin non-empty

References (host + all three viruses) are emitted by
``write_references.py`` so the same sequences feed bwa, kraken2-build,
kaiju-mkbwt, makeblastdb and the virus parquet.

The output FASTQ files are gzipped when the requested path ends in
``.gz``.
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

QUAL = "I"  # Phred 40, ASCII 73


def _pseudo_random_seq(seed: int, length: int) -> str:
    """Deterministic pseudo-random DNA of `length` from `seed`. Used so the
    smoke fixtures are reproducible without committing 10 kb of FASTA.
    """
    rng = random.Random(seed)
    bases = "ACGT"
    return "".join(rng.choice(bases) for _ in range(length))


# Single 5 kb host. The same host is used for every sample so the
# host-removal stage exercises consistently across the batch.
HOST_REF = _pseudo_random_seq(seed=1, length=5000)


# Three distinct viral references. Names and taxids are stable so the
# kraken/kaiju/blast/parquet builders can hard-code them. Sequences are
# fully pseudo-random (non-repetitive) so MEGAHIT has a path through the
# de Bruijn graph; each seed is unique so the three viruses do not share
# k-mers and the classifiers can tell them apart.
VIRUS_REFS: dict[str, dict] = {
    "alpha": {"taxid": 100001, "sequence": _pseudo_random_seq(seed=2, length=5000)},
    "beta": {"taxid": 100002, "sequence": _pseudo_random_seq(seed=3, length=5000)},
    "gamma": {"taxid": 100003, "sequence": _pseudo_random_seq(seed=4, length=5000)},
}


def make_read(template: str, start: int, length: int = 150) -> str:
    return template[start : start + length].ljust(length, "N")


def write_pair(
    out_r1: Path,
    out_r2: Path,
    *,
    sample: str,
    virus: str,
    seed: int,
    n_host: int = 60,
    n_virus: int = 800,
    n_random: int = 50,
) -> None:
    if virus not in VIRUS_REFS:
        raise SystemExit(
            f"Unknown virus '{virus}'. Known: {sorted(VIRUS_REFS.keys())}"
        )

    rng = random.Random(seed)
    virus_seq = VIRUS_REFS[virus]["sequence"]

    records: list[tuple[str, str, str]] = []  # (id, r1, r2)
    for i in range(n_host):
        pos = rng.randint(0, len(HOST_REF) - 400)
        r1 = make_read(HOST_REF, pos)
        r2 = make_read(HOST_REF, pos + 200)
        records.append((f"{sample}_host_{i:03d}", r1, r2))
    for i in range(n_virus):
        pos = rng.randint(0, len(virus_seq) - 400)
        r1 = make_read(virus_seq, pos)
        r2 = make_read(virus_seq, pos + 200)
        records.append((f"{sample}_{virus}_{i:03d}", r1, r2))
    bases = "ACGT"
    for i in range(n_random):
        r1 = "".join(rng.choice(bases) for _ in range(150))
        r2 = "".join(rng.choice(bases) for _ in range(150))
        records.append((f"{sample}_random_{i:03d}", r1, r2))

    rng.shuffle(records)

    opener = gzip.open if out_r1.suffix == ".gz" else open
    with opener(out_r1, "wt") as f1, opener(out_r2, "wt") as f2:
        for read_id, r1, r2 in records:
            f1.write(f"@{read_id}/1\n{r1}\n+\n{QUAL * len(r1)}\n")
            f2.write(f"@{read_id}/2\n{r2}\n+\n{QUAL * len(r2)}\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-r1", type=Path, required=True)
    p.add_argument("--out-r2", type=Path, required=True)
    p.add_argument("--sample", default="sample")
    p.add_argument(
        "--virus",
        required=True,
        choices=sorted(VIRUS_REFS.keys()),
        help="Which of the three synthetic viruses this sample carries.",
    )
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--n-virus", type=int, default=800)
    p.add_argument("--n-host", type=int, default=60)
    p.add_argument("--n-random", type=int, default=50)
    args = p.parse_args()

    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r2.parent.mkdir(parents=True, exist_ok=True)
    write_pair(
        args.out_r1,
        args.out_r2,
        sample=args.sample,
        virus=args.virus,
        seed=args.seed,
        n_host=args.n_host,
        n_virus=args.n_virus,
        n_random=args.n_random,
    )


def get_host_reference() -> str:
    return HOST_REF


def get_virus_references() -> dict[str, dict]:
    """Return the dict of {name: {taxid, sequence}} so the reference and
    DB-builder scripts share one source of truth.
    """
    return VIRUS_REFS


if __name__ == "__main__":
    main()
