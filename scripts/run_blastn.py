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
from subprocess import CalledProcessError
import sys

# When invoked as ``script:`` Snakemake puts the rule's working
# directory on sys.path but not the project root. Add the script
# directory's parent so ``from scripts.functions`` resolves.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.functions import run_blastn  # noqa: E402

snakemake = snakemake  # type: ignore[name-defined]  # injected at runtime

Path(snakemake.output.blast).parent.mkdir(parents=True, exist_ok=True)
try:
    df = run_blastn(
        contigs_csv=snakemake.input.contigs,
        db=snakemake.params.db,
        temp_file=snakemake.params.temp_file,
        threads=snakemake.threads,
    )
except CalledProcessError as exc:
    # run_blastn captures blastn's stderr on the exception, but it is
    # otherwise invisible (the rule's log stays empty, which makes a DB /
    # query failure very hard to diagnose). Write the real error to the
    # rule log and stderr before failing.
    cmd = " ".join(exc.cmd) if isinstance(exc.cmd, list) else str(exc.cmd)
    detail = (
        f"blastn failed (exit {exc.returncode})\n"
        f"command: {cmd}\n"
        f"--- blastn stderr ---\n{exc.stderr or ''}\n"
        f"--- blastn stdout ---\n{exc.stdout or ''}\n"
    )
    for logpath in getattr(snakemake, "log", []) or []:
        try:
            with open(logpath, "a") as fh:
                fh.write(detail)
        except OSError:
            pass
    sys.stderr.write(detail)
    raise
df.to_csv(snakemake.output.blast, index=False)

temp = Path(snakemake.params.temp_file)
if temp.exists():
    try:
        os.remove(temp)
    except OSError as exc:
        logging.warning("Could not remove BLASTN temp file %s: %s", temp, exc)
