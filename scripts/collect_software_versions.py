"""Merge per-env conda-meta dumps into a resolved software-version table.

Under ``--sdm conda`` each tool rule runs inside a materialised conda
prefix. A thin probe rule that references the same ``envs/<env>.yaml``
reuses that exact prefix and dumps the ``conda-meta`` json filename
stems (``<name>-<version>-<build>``) to ``logs/versions/<env>.tsv``, one
``<env>\\t<stem>`` line per package. This module parses those stems and
merges them into a single ``software_versions.tsv`` plus a headline
``{package: version}`` map -- the resolved versions of exactly what ran,
which the provenance sidecar, the run-info CSV and the HTML report reuse.

Standard library only (parsing must not assume pandas): a probe env such
as ``fastp`` carries no Python, so parsing lives here, not in the shell.

Note: no ``from __future__ import annotations`` here. Snakemake prepends
its own boilerplate when materialising a ``script:`` rule, so a future
import would no longer sit at the top of the file and would raise a
SyntaxError. The py>=3.11 rule envs evaluate the PEP 604/585 annotations
below natively.
"""

import csv
from pathlib import Path

# One headline tool per probed env. The flat table keeps every package;
# the headline map surfaces the tool the env exists for (e.g. the
# ``blastn`` env's headline package is ``blast``, the ``kraken`` env's is
# ``kraken2``), skipping incidental dependencies like python or libgcc.
ENV_PRIMARY_PACKAGE: dict[str, str] = {
    "fastp": "fastp",
    "bwa": "bwa",
    "samtools": "samtools",
    "kraken": "kraken2",
    "kaiju": "kaiju",
    "megahit": "megahit",
    "spades": "spades",
    "pilon": "pilon",
    "blastn": "blast",
    "checkv": "checkv",
    "mosdepth": "mosdepth",
    "multiqc": "multiqc",
    "genomad": "genomad",
    "quast": "quast",
    "hostile": "hostile",
}


def parse_conda_meta_stem(stem: str) -> tuple[str, str, str]:
    """Split a ``conda-meta`` filename stem into ``(name, version, build)``.

    The stem is ``<name>-<version>-<build>``; the build and version are
    the last two hyphen-delimited fields, so the remainder -- which may
    itself contain hyphens, e.g. ``perl-archive-tar`` -- is the package
    name. A stem without two hyphens cannot be split, so it is returned
    whole as the name with empty version/build rather than dropped.
    """
    rest, sep, build = stem.rpartition("-")
    if not sep:
        return stem, "", ""
    name, sep, version = rest.rpartition("-")
    if not sep:
        return stem, "", ""
    return name, version, build


def collect_versions(probe_tsvs: list[Path]) -> list[dict[str, str]]:
    """Read every per-env probe TSV and return flat version rows.

    Each row is ``{env, package, version, build}``. Missing probe files
    are skipped so an optional stage that did not run leaves no rows
    rather than failing the merge.
    """
    rows: list[dict[str, str]] = []
    for tsv in probe_tsvs:
        path = Path(tsv)
        if not path.is_file():
            continue
        with path.open() as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                env, _, stem = line.partition("\t")
                if not stem:
                    continue
                name, version, build = parse_conda_meta_stem(stem)
                rows.append(
                    {"env": env, "package": name, "version": version, "build": build}
                )
    rows.sort(key=lambda r: (r["env"], r["package"]))
    return rows


def read_versions_table(path: Path) -> list[dict[str, str]]:
    """Read a merged ``software_versions.tsv`` (header + env/package/
    version/build rows) back into row dicts.

    This is the inverse of :func:`write_software_versions`; use it to
    reload the collector's output, NOT :func:`collect_versions` (which
    parses the raw ``env<TAB>stem`` probe dumps). Returns an empty list
    when the file is absent so a legacy run degrades cleanly.
    """
    p = Path(path)
    if not p.is_file():
        return []
    rows: list[dict[str, str]] = []
    with p.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            if len(values) == len(header):
                rows.append(dict(zip(header, values)))
    return rows


def headline_versions(rows: list[dict[str, str]]) -> dict[str, str]:
    """Reduce flat rows to a ``{package: version}`` map of headline tools.

    For each env only its ``ENV_PRIMARY_PACKAGE`` entry is kept, keyed by
    the package name (e.g. ``{"kraken2": "2.1.3"}``), so the map reads as
    the versions of the tools the pipeline actually invoked.
    """
    out: dict[str, str] = {}
    for row in rows:
        primary = ENV_PRIMARY_PACKAGE.get(row["env"])
        if primary and row["package"] == primary:
            out[primary] = row["version"]
    return out


def write_software_versions(rows: list[dict[str, str]], out_path: Path) -> None:
    """Write the flat version rows as a tab-separated file with a header."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["env", "package", "version", "build"])
        for row in rows:
            writer.writerow([row["env"], row["package"], row["version"], row["build"]])


def main() -> None:
    # snakemake is injected into globals by the Snakemake `script:` runner.
    sm = globals()["snakemake"]  # noqa: F821 (provided at runtime)
    probe_tsvs = [Path(p) for p in sm.input]
    rows = collect_versions(probe_tsvs)
    write_software_versions(rows, Path(sm.output[0]))


if __name__ == "__main__":
    main()
