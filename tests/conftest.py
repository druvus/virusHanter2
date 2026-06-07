"""Shared pytest fixtures for the virusHanter2 test suite."""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _isolate_taxdump_cache(tmp_path_factory, monkeypatch):
    """Redirect the taxdump pickle cache to a throwaway dir for each test.

    ``scripts/functions`` caches parsed nodes.dmp / names.dmp structures
    under ``$XDG_CACHE_HOME`` (falling back to the OS temp dir). Pinning
    it to a fresh per-test directory keeps the suite hermetic: tests that
    parse taxdump fixtures never collide on a shared cache and never
    leave cache files in the developer's real cache dir.
    """
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path_factory.mktemp("xdg_cache")))
