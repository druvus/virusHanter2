"""Lightweight, dependency-free helpers shared by the assembler scripts.

These run inside the minimal per-rule conda envs (``envs/megahit.yaml``,
``envs/spades.yaml``) which carry only the assembler binary and no
scientific Python stack. They must therefore import nothing beyond the
standard library; in particular they must not pull in ``scripts.functions``,
which imports pandas and numpy and would raise ModuleNotFoundError in those
envs. ``scripts.functions`` re-exports both helpers for driver-env callers
and the test-suite.
"""

from pathlib import Path


def assembler_max_attempts(cfg: dict, is_target_platform: bool) -> int:
    """Return the total number of assembler attempts to make.

    On the target platform, reads ``ASSEMBLER_RETRIES`` from the config
    (or the legacy ``MEGAHIT_RETRIES`` key for backward compatibility)
    and returns that value plus one for the initial attempt.  On all
    other platforms a single attempt is made.

    Using a generic key means both MEGAHIT and SPAdes share the same
    retry budget on Apple Silicon, where non-deterministic SIGSEGV or
    library-size failures can occur.
    """
    if not is_target_platform:
        return 1
    n = int(cfg.get("ASSEMBLER_RETRIES", cfg.get("MEGAHIT_RETRIES", 4)))
    return n + 1


def write_dummy_contig(path: str) -> None:
    """Write a minimal placeholder FASTA contig to ``path``.

    Used when an assembler produces no usable output so that downstream
    rules always receive a syntactically valid FASTA.  The contig
    carries a synthetic 200 bp sequence and will not match any real
    database entry.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(">DUMMY_CONTIG\n")
        fh.write("TTAACCTTGG" * 20 + "\n")
