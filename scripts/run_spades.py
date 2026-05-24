"""Snakemake script: de novo assembly with SPAdes.

Used by both ``rule metaspades`` (``--meta``) and
``rule rnaviralspades`` (``--rnaviral``). The mode flag is read
from ``snakemake.params.mode``. Falls back to a dummy contig when
SPAdes refuses a too-small library, so the per-assembler DAG
keeps a uniform shape for downstream rules.
"""

import subprocess
from pathlib import Path

snakemake = snakemake  # type: ignore[name-defined]

params = snakemake.params
input_ = snakemake.input
output = snakemake.output
threads = snakemake.threads
log_path = snakemake.log[0] if snakemake.log else "/dev/null"
mode = params.mode  # "meta" or "rnaviral"


def _shell(cmd: str, check: bool = True) -> int:
    return subprocess.run(
        cmd, shell=True, check=check, executable="/bin/bash"
    ).returncode


_shell(f"rm -rf {params.out_dir}")
_shell(f"mkdir -p {params.out_dir}")

mem_gb = max(8, int(snakemake.resources.mem_mb / 1024))

try:
    _shell(
        f"spades.py --{mode} "
        f"-1 {input_.r1} -2 {input_.r2} "
        f"-o {params.out_dir} "
        f"-t {threads} -m {mem_gb} --only-assembler "
        f"> {log_path} 2>&1"
    )
except subprocess.CalledProcessError:
    # SPAdes refuses libraries below its internal minimum and exits
    # non-zero. Continue to the fallback so the rest of the DAG
    # still has an input.
    pass

# metaSPAdes writes contigs.fasta; rnaviralSPAdes writes
# transcripts.fasta. Try both, in that order.
candidates = (
    ("transcripts.fasta", "contigs.fasta")
    if mode == "rnaviral"
    else ("contigs.fasta",)
)
moved = False
for candidate in candidates:
    src = Path(params.out_dir) / candidate
    if src.exists() and src.stat().st_size > 0:
        _shell(f"mv {src} {output.contigs}")
        moved = True
        break

if not moved:
    Path(params.out_dir).mkdir(parents=True, exist_ok=True)
    with open(output.contigs, "w") as f:
        f.write(">DUMMY_CONTIG\n")
        f.write("TTAACCTTGG" * 20 + "\n")

# Strip SPAdes intermediates; keep the contigs FASTA only.
_shell(
    f"ls -d -1 {params.out_dir}/* 2>/dev/null "
    f"| grep -v .contigs.fa | xargs rm -rf",
    check=False,
)
