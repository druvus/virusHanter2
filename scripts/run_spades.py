"""Snakemake script: de novo assembly with SPAdes.

Used by both ``rule metaspades`` (``--meta``) and
``rule rnaviralspades`` (``--rnaviral``). The mode flag is read
from ``snakemake.params.mode``. Retries on Apple Silicon (where
non-deterministic library-size failures occur) using the same
configurable retry budget as MEGAHIT (``ASSEMBLER_RETRIES`` /
legacy ``MEGAHIT_RETRIES``). Falls back to a dummy contig when
all attempts fail so the per-assembler DAG keeps a uniform shape.
"""

import platform
import shlex
import shutil
import subprocess
from pathlib import Path

from scripts.assembler_utils import assembler_max_attempts, write_dummy_contig

snakemake = snakemake  # type: ignore[name-defined]

params = snakemake.params
input_ = snakemake.input
output = snakemake.output
threads = snakemake.threads
log_path = snakemake.log[0] if snakemake.log else "/dev/null"
mode = params.mode  # "meta" or "rnaviral"


# Resolve bash at import time so non-standard installations are handled
# transparently. Falls back to /bin/bash if bash is not found on PATH.
_BASH = shutil.which("bash") or "/bin/bash"


def _shell(cmd: str, check: bool = True) -> int:
    return subprocess.run(
        cmd, shell=True, check=check, executable=_BASH
    ).returncode


is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
max_attempts = assembler_max_attempts(snakemake.config, is_apple_silicon)

mem_gb = max(8, int(snakemake.resources.mem_mb / 1024))

# metaSPAdes writes contigs.fasta; rnaviralSPAdes writes
# transcripts.fasta. Try both, in that order.
candidates = (
    ("transcripts.fasta", "contigs.fasta")
    if mode == "rnaviral"
    else ("contigs.fasta",)
)

# Shell-quote every interpolated path: a RESULTS_FOLDER on a macOS
# external volume (e.g. "/Volumes/My Passport/...") contains spaces that
# would otherwise split into multiple shell arguments.
_q_out_dir = shlex.quote(str(params.out_dir))
_q_r1 = shlex.quote(str(input_.r1))
_q_r2 = shlex.quote(str(input_.r2))
_q_contigs = shlex.quote(str(output.contigs))
_q_log = shlex.quote(str(log_path))

success = False
for attempt in range(1, max_attempts + 1):
    _shell(f"rm -rf {_q_out_dir}", check=False)
    _shell(f"mkdir -p {_q_out_dir}")

    try:
        _shell(
            f"spades.py --{mode} "
            f"-1 {_q_r1} -2 {_q_r2} "
            f"-o {_q_out_dir} "
            f"-t {threads} -m {mem_gb} --only-assembler "
            f"> {_q_log} 2>&1"
        )
    except subprocess.CalledProcessError:
        # SPAdes refuses libraries below its internal minimum and exits
        # non-zero. Continue to the next attempt or fall back to the
        # dummy contig.
        continue

    for candidate in candidates:
        src = Path(params.out_dir) / candidate
        if src.exists() and src.stat().st_size > 0:
            _shell(f"mv {shlex.quote(str(src))} {_q_contigs}")
            success = True
            break

    if success:
        break

if not success:
    write_dummy_contig(output.contigs)

# Strip SPAdes intermediates; keep the contigs FASTA(s) only. Done in
# Python rather than `ls | grep -v .contigs.fa | xargs rm -rf`: the old
# pipe used an unescaped, unanchored pattern and an `xargs` without `-r`,
# so an empty pipe could run `rm -rf` with no operand on some platforms.
# The contigs were renamed to `{sample}.contigs.fa` (suffix `.fa`) above;
# every other SPAdes output (`*.fasta`, `*.gfa`, K-mer dirs) is dropped.
out_dir = Path(params.out_dir)
if out_dir.is_dir():
    for entry in out_dir.iterdir():
        if entry.suffix == ".fa":
            continue
        if entry.is_dir():
            shutil.rmtree(entry, ignore_errors=True)
        else:
            entry.unlink(missing_ok=True)
