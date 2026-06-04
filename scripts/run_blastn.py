"""Snakemake script: run BLASTN on contigs.

Receives the standard ``snakemake`` magic object. Delegates the
actual work to ``run_blastn`` in ``scripts/functions.py`` so the
existing implementation (and its tests) keeps a single source of
truth. Runs in the rule's conda env (``envs/blastn.yaml``) so the
``blastn`` binary is on PATH.
"""

import logging
import os
from pathlib import Path
import sys

# When invoked as ``script:`` Snakemake puts the rule's working
# directory on sys.path but not the project root. Add the script
# directory's parent so ``from scripts.functions`` resolves.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.functions import run_blastn  # noqa: E402

snakemake = snakemake  # type: ignore[name-defined]  # injected at runtime

Path(snakemake.output.blast).parent.mkdir(parents=True, exist_ok=True)
df = run_blastn(
    contigs_csv=snakemake.input.contigs,
    db=snakemake.params.db,
    temp_file=snakemake.params.temp_file,
    threads=snakemake.threads,
)
df.to_csv(snakemake.output.blast, index=False)

temp = Path(snakemake.params.temp_file)
if temp.exists():
    try:
        os.remove(temp)
    except OSError as exc:
        logging.warning("Could not remove BLASTN temp file %s: %s", temp, exc)
