"""Snakemake script: de novo assembly with MEGAHIT.

Runs in the rule's conda env (``envs/megahit.yaml``) so the
``megahit`` binary is on PATH. Carries the Apple Silicon
work-arounds (``--no-hw-accel``, thread cap, k-range clamp, retry
loop) and writes a dummy contig when all attempts fail so the
per-assembler DAG keeps a uniform shape.
"""

import platform
import shutil
import subprocess
from pathlib import Path

from scripts.functions import assembler_max_attempts, write_dummy_contig

snakemake = snakemake  # type: ignore[name-defined]

config = snakemake.config
params = snakemake.params
input_ = snakemake.input
output = snakemake.output
threads = snakemake.threads
log_path = snakemake.log[0] if snakemake.log else "/dev/null"
sample = snakemake.wildcards.sample


# Resolve the bash executable at import time. ``shutil.which`` searches PATH,
# which covers non-standard installations (e.g. /usr/local/bin/bash on macOS
# with Homebrew). Fall back to /bin/bash if bash is not on PATH; the subprocess
# call will then fail loudly if that path also does not exist.
_BASH = shutil.which("bash") or "/bin/bash"


def _shell(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, executable=_BASH)


def _shell_tolerate(cmd: str) -> bool:
    """Run a command, return True on success and False on non-zero exit."""
    return (
        subprocess.run(cmd, shell=True, check=False, executable=_BASH).returncode
        == 0
    )


_shell(f"rm -rf {params.out_dir}")

mem_fraction = float(config.get("MEGAHIT_MEM_FRACTION", 0.5))
is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
no_hw_accel = "--no-hw-accel " if is_apple_silicon else ""
# On Apple Silicon (Darwin/arm64) MEGAHIT SIGSEGVs with the default
# k-min of 21 when the input library is small (a known upstream bug).
# k-min 27 is the smallest value that avoids the crash on tested inputs.
# k-max 57 is chosen as a stable upper bound: values above 63 require
# the MEGAHIT large-k binary (not shipped by bioconda), and 57 is
# sufficient for typical viral contig recovery from short reads.
kmin_flag = "--k-min 27 " if is_apple_silicon else ""
kmax_flag = "--k-max 57 " if is_apple_silicon else ""
mh_threads = min(threads, 2) if is_apple_silicon else threads

max_attempts = assembler_max_attempts(config, is_apple_silicon)

success = False
for attempt in range(1, max_attempts + 1):
    _shell(f"rm -rf {params.out_dir}")
    try:
        _shell(
            "megahit "
            f"-1 {input_.r1} -2 {input_.r2} "
            f"-o {params.out_dir} "
            f"--out-prefix {sample} "
            f"-t {mh_threads} "
            f"-m {mem_fraction} "
            f"{no_hw_accel}{kmin_flag}{kmax_flag}"
            f"2> {log_path}"
        )
        if Path(output.contigs).exists() and Path(output.contigs).stat().st_size > 0:
            success = True
            break
    except subprocess.CalledProcessError:
        continue

if not success:
    write_dummy_contig(output.contigs)

# Drop MEGAHIT intermediates; keep the contigs FASTA only.
_shell_tolerate(
    f"ls -d -1 {params.out_dir}/* 2>/dev/null | grep -v .fa | xargs rm -rf"
)

if Path(output.contigs).read_text() == "":
    write_dummy_contig(output.contigs)
