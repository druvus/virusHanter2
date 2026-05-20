"""
Emit the viral-reference Parquet (`name`, `sequence`, `tax_id`) used by
the `bwa_align_to_kraken_hits` rule. Mirrors the three synthetic viruses
defined in synthesize_fastq.py.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from synthesize_fastq import get_virus_references


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    rows = [
        {"name": name, "sequence": info["sequence"], "tax_id": info["taxid"]}
        for name, info in get_virus_references().items()
    ]
    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)


if __name__ == "__main__":
    main()
