"""
Emit the viral-reference Parquet (`name`, `sequence`, `tax_id`) used by the
`bwa_align_to_kraken_hits` rule.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from synthesize_fastq import get_virus_reference


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    df = pd.DataFrame(
        [
            {
                "name": "synthetic_virus",
                "sequence": get_virus_reference(),
                "tax_id": 100001,
            },
            {
                "name": "synthetic_virus_b",
                "sequence": get_virus_reference()[::-1],
                "tax_id": 100002,
            },
        ]
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)


if __name__ == "__main__":
    main()
