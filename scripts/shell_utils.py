"""Dependency-free shell helpers shared by the pipeline's script: rules.

Standard library only, so the helpers import cleanly from any per-rule
conda env (mirroring ``assembler_utils``). They must not pull in
``scripts.functions`` or any scientific-Python dependency.
"""

from __future__ import annotations

import shutil
import subprocess


def resolve_bash() -> str:
    """Return the bash executable path, falling back to ``/bin/bash``.

    ``shutil.which`` covers non-standard installations (e.g. Homebrew's
    ``/usr/local/bin/bash`` on macOS); the fallback lets the subprocess
    call fail loudly if that path is also absent.
    """
    return shutil.which("bash") or "/bin/bash"


def run_piped(cmd: str, *, bash: str | None = None, check: bool = True) -> int:
    """Run ``cmd`` under bash with ``pipefail`` enabled; return the exit code.

    Without ``pipefail`` a shell pipeline reports only the last stage's
    exit status, so a mid-pipe failure -- for example ``bwa mem``
    aborting (out of memory, truncated FASTQ, index mismatch) while the
    downstream ``samtools sort`` still exits 0 -- goes undetected and a
    truncated-but-valid output file looks like a clean success. Enabling
    ``pipefail`` makes the pipeline exit non-zero when any stage fails,
    so ``check=True`` surfaces the failure instead of letting it pass
    silently downstream.
    """
    executable = bash or resolve_bash()
    return subprocess.run(
        f"set -o pipefail; {cmd}",
        shell=True,
        check=check,
        executable=executable,
    ).returncode
