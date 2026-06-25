"""Tests for the shared pipefail shell helper.

A bash pipeline without ``pipefail`` reports only the last stage's exit
status, so a mid-pipe failure (e.g. ``bwa mem`` aborting while
``samtools sort`` still exits 0) goes undetected and a truncated BAM
looks like success. ``run_piped`` enables ``pipefail`` so ``check=True``
surfaces the failure.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.shell_utils import run_piped  # noqa: E402


def test_run_piped_detects_first_stage_failure():
    # `false | true`: the last stage (true) exits 0, so without pipefail
    # the pipeline would report success. With pipefail the failing first
    # stage propagates and check=True raises.
    with pytest.raises(subprocess.CalledProcessError):
        run_piped("false | true")


def test_run_piped_succeeds_when_all_stages_pass():
    assert run_piped("true | true") == 0


def test_run_piped_no_check_returns_nonzero_without_raising():
    rc = run_piped("false | true", check=False)
    assert rc != 0
