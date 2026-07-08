#!/usr/bin/env python3
"""Export the per-sample HTML reports from a virusHanter2 results tree
into a separate output location, preserving the per-run folder layout.

A virusHanter2 results root holds one folder per sequencing run::

    <results>/<run>/<sample>/REPORT/<sample>.html
    <results>/<run>/run_information_<run>.csv
    <results>/<run>/per_virus_<run>.csv
    <results>/<run>/run_provenance_<run>.json
    <results>/<run>/run_provenance_<run>.tsv
    <results>/<run>/software_versions.tsv

The exporter mirrors each run folder into the destination, copying every
sample HTML report and -- unless ``--reports-only`` is given -- the small
run-level summary files that accompany them. By default the reports are
flattened to ``<dest>/<run>/<sample>.html`` (the intermediate
``<sample>/REPORT/`` directories carry no information for a share-out);
pass ``--preserve-tree`` to reproduce the full source path instead.

Only the copied files are touched; the source tree is never modified.

Examples
--------

Export every analysed run under the current results root::

    python scripts/export_reports.py ~/vh2-results ~/vh2-reports-export

Preview without copying, HTML reports only::

    python scripts/export_reports.py ~/vh2-results ~/out --reports-only --dry-run

Export two named runs, keeping the original folder tree::

    python scripts/export_reports.py ~/vh2-results ~/out --preserve-tree \\
        --run 260611_M00568_0769_000000000-DRMJL \\
        --run 251015_M00568_0723_000000000-DRRKK
"""
from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Small run-level summary files that live directly in each run folder and
# usefully travel with the reports. Globbed per run; any that are absent
# are simply skipped.
SUMMARY_GLOBS = (
    "run_information_*.csv",
    "per_virus_*.csv",
    "run_provenance_*.json",
    "run_provenance_*.tsv",
    "software_versions.tsv",
)


@dataclass
class RunExport:
    """The files discovered for one run, ready to copy."""

    name: str
    reports: list[Path] = field(default_factory=list)
    summaries: list[Path] = field(default_factory=list)


def find_sample_reports(run_dir: Path) -> list[Path]:
    """Return the per-sample HTML reports under a run folder.

    Matches ``<sample>/REPORT/<sample>.html``. Dot-prefixed entries
    (macOS AppleDouble '._*' shadows and other hidden files) are ignored
    at every level.
    """
    reports: list[Path] = []
    for sample_dir in sorted(run_dir.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name.startswith("."):
            continue
        report_dir = sample_dir / "REPORT"
        if not report_dir.is_dir():
            continue
        reports.extend(
            sorted(
                html
                for html in report_dir.glob("*.html")
                if not html.name.startswith(".")
            )
        )
    return reports


def find_summaries(run_dir: Path) -> list[Path]:
    """Return the run-level summary files present in a run folder."""
    found: list[Path] = []
    for pattern in SUMMARY_GLOBS:
        found.extend(
            sorted(p for p in run_dir.glob(pattern) if not p.name.startswith("."))
        )
    return found


def discover_runs(source: Path, only: set[str] | None) -> list[RunExport]:
    """Scan the results root for analysed runs.

    A run is any immediate subdirectory carrying at least one sample HTML
    report. When ``only`` is given, restrict discovery to those run
    names.
    """
    runs: list[RunExport] = []
    for run_dir in sorted(source.iterdir()):
        if not run_dir.is_dir() or run_dir.name.startswith("."):
            continue
        if only is not None and run_dir.name not in only:
            continue
        reports = find_sample_reports(run_dir)
        if not reports:
            continue
        runs.append(
            RunExport(
                name=run_dir.name,
                reports=reports,
                summaries=find_summaries(run_dir),
            )
        )
    return runs


def destination_for(
    source_file: Path, run_dir: Path, dest_run: Path, preserve_tree: bool
) -> Path:
    """Work out where a source file lands in the destination run folder."""
    if preserve_tree:
        return dest_run / source_file.relative_to(run_dir)
    return dest_run / source_file.name


def copy_file(src: Path, dst: Path, *, skip_existing: bool, dry_run: bool) -> str:
    """Copy one file, returning a one-word status for the summary tally."""
    if skip_existing and dst.exists():
        return "skipped"
    if dry_run:
        return "planned"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return "copied"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("source", type=Path, help="virusHanter2 results root to read from")
    p.add_argument("dest", type=Path, help="output folder to copy reports into")
    p.add_argument(
        "--run",
        dest="runs",
        action="append",
        metavar="RUN_NAME",
        help="export only this run folder; repeatable. Default: all analysed runs.",
    )
    p.add_argument(
        "--reports-only",
        action="store_true",
        help="copy only the HTML reports, not the run-level summary files.",
    )
    p.add_argument(
        "--preserve-tree",
        action="store_true",
        help="keep the <sample>/REPORT/ path instead of flattening to <sample>.html.",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="do not overwrite a file that already exists in the destination.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="report what would be copied without writing anything.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if not args.source.is_dir():
        print(f"error: source is not a directory: {args.source}", file=sys.stderr)
        return 2

    only = set(args.runs) if args.runs else None
    runs = discover_runs(args.source, only)

    if only:
        missing = only - {r.name for r in runs}
        for name in sorted(missing):
            print(f"warning: no reports found for run '{name}'", file=sys.stderr)

    if not runs:
        print("No analysed runs with reports were found; nothing to export.")
        return 1

    tally = {"copied": 0, "skipped": 0, "planned": 0}
    for run in runs:
        run_dir = args.source / run.name
        dest_run = args.dest / run.name
        files = list(run.reports)
        if not args.reports_only:
            files += run.summaries
        for src in files:
            dst = destination_for(src, run_dir, dest_run, args.preserve_tree)
            status = copy_file(
                src, dst, skip_existing=args.skip_existing, dry_run=args.dry_run
            )
            tally[status] += 1
        print(
            f"{run.name}: {len(run.reports)} report(s)"
            + ("" if args.reports_only else f", {len(run.summaries)} summary file(s)")
        )

    verb = "Would copy" if args.dry_run else "Copied"
    line = f"{verb} {tally['copied'] or tally['planned']} file(s) from {len(runs)} run(s)"
    if tally["skipped"]:
        line += f"; skipped {tally['skipped']} already present"
    print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
