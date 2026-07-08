"""Tests for scripts/export_reports.py.

The exporter mirrors a virusHanter2 results tree into a destination,
copying the per-sample HTML reports (flattened by default) and the
run-level summary files, without touching the source.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.export_reports import discover_runs, main  # noqa: E402


def _make_run(root: Path, run: str, samples: list[str]) -> Path:
    run_dir = root / run
    for s in samples:
        rep = run_dir / s / "REPORT"
        rep.mkdir(parents=True)
        (rep / f"{s}.html").write_text(f"<html>{s}</html>")
    (run_dir / f"run_information_{run}.csv").write_text("col\n1\n")
    (run_dir / f"per_virus_{run}.csv").write_text("virus\n")
    (run_dir / f"run_provenance_{run}.json").write_text("{}")
    (run_dir / "software_versions.tsv").write_text("tool\tversion\n")
    return run_dir


def test_discover_and_export_flattened(tmp_path: Path) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R", "s2_R"])
    dest = tmp_path / "out"

    rc = main([str(src), str(dest)])

    assert rc == 0
    assert (dest / "RUN_A" / "s1_R.html").read_text() == "<html>s1_R</html>"
    assert (dest / "RUN_A" / "s2_R.html").exists()
    # Summaries travel by default.
    assert (dest / "RUN_A" / "run_information_RUN_A.csv").exists()
    assert (dest / "RUN_A" / "per_virus_RUN_A.csv").exists()
    assert (dest / "RUN_A" / "software_versions.tsv").exists()
    # Source is untouched.
    assert (src / "RUN_A" / "s1_R" / "REPORT" / "s1_R.html").exists()


def test_preserve_tree(tmp_path: Path) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R"])
    dest = tmp_path / "out"

    main([str(src), str(dest), "--preserve-tree"])

    assert (dest / "RUN_A" / "s1_R" / "REPORT" / "s1_R.html").exists()
    assert not (dest / "RUN_A" / "s1_R.html").exists()


def test_reports_only_skips_summaries(tmp_path: Path) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R"])
    dest = tmp_path / "out"

    main([str(src), str(dest), "--reports-only"])

    assert (dest / "RUN_A" / "s1_R.html").exists()
    assert not (dest / "RUN_A" / "run_information_RUN_A.csv").exists()


def test_ignores_appledouble_and_sampleless_dirs(tmp_path: Path) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R"])
    # AppleDouble shadow next to a real report; a failed sample with no
    # REPORT dir; a hidden run folder. None should be exported.
    (src / "RUN_A" / "s1_R" / "REPORT" / "._s1_R.html").write_text("junk")
    (src / "RUN_A" / "s2_failed").mkdir()
    _make_run(src, ".hidden_run", ["x_R"])
    dest = tmp_path / "out"

    runs = discover_runs(src, None)
    assert [r.name for r in runs] == ["RUN_A"]
    assert len(runs[0].reports) == 1

    main([str(src), str(dest)])
    assert not (dest / "RUN_A" / "._s1_R.html").exists()
    assert not (dest / ".hidden_run").exists()


def test_run_filter_and_missing_warning(tmp_path: Path, capsys) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R"])
    _make_run(src, "RUN_B", ["s1_R"])
    dest = tmp_path / "out"

    main([str(src), str(dest), "--run", "RUN_A", "--run", "RUN_MISSING"])

    assert (dest / "RUN_A" / "s1_R.html").exists()
    assert not (dest / "RUN_B").exists()
    err = capsys.readouterr().err
    assert "RUN_MISSING" in err


def test_dry_run_writes_nothing(tmp_path: Path) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R"])
    dest = tmp_path / "out"

    rc = main([str(src), str(dest), "--dry-run"])

    assert rc == 0
    assert not dest.exists()


def test_skip_existing(tmp_path: Path) -> None:
    src = tmp_path / "results"
    _make_run(src, "RUN_A", ["s1_R"])
    dest = tmp_path / "out"
    main([str(src), str(dest)])

    # Mutate the destination, re-export with --skip-existing: unchanged.
    target = dest / "RUN_A" / "s1_R.html"
    target.write_text("EDITED")
    main([str(src), str(dest), "--skip-existing"])
    assert target.read_text() == "EDITED"

    # Without the flag, it is overwritten from source.
    main([str(src), str(dest)])
    assert target.read_text() == "<html>s1_R</html>"
