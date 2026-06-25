"""Tests for the taxdump parse cache in ``scripts/functions.py``.

``parse_nodes_dmp`` and ``_load_taxdump_for_species_walkup`` pickle
their parsed structures to a per-user cache dir keyed by source path +
mtime + size, because the full NCBI taxdump (hundreds of MB) would
otherwise be re-parsed once per rule invocation. These tests lock in the
contract:

  - a cache hit returns a result identical to a cold parse,
  - the cache is actually consulted (a poisoned cache is returned),
  - a changed source file (new mtime/size) invalidates the cache,
  - a corrupt cache file degrades to a fresh parse rather than raising,
  - a missing source file is not cached,
  - the cache key separates distinct sources and tracks file size,
  - the low-level store/load helpers round-trip and fail silently.

The cache directory is redirected to an isolated ``tmp_path`` via
``XDG_CACHE_HOME`` so the tests never touch the real cache.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import functions as F  # noqa: E402

# Minimal NCBI .dmp fixtures. nodes.dmp fields are
# ``tax_id | parent | rank | ...`` joined by ``\t|\t`` and terminated by
# ``\t|``; names.dmp fields are ``tax_id | name | unique | name_class``.
NODES = (
    "1\t|\t1\t|\tno rank\t|\n"
    "10239\t|\t1\t|\tacellular root\t|\n"
    "2559587\t|\t10239\t|\trealm\t|\n"
    "3050294\t|\t2559587\t|\tgenus\t|\n"
    "3050299\t|\t3050294\t|\tspecies\t|\n"
    "10376\t|\t3050299\t|\tno rank\t|\n"
)

NAMES = (
    "3050299\t|\tLymphocryptovirus humangamma4\t|\t\t|\tscientific name\t|\n"
    "10376\t|\thuman gammaherpesvirus 4\t|\t\t|\tscientific name\t|\n"
    "10376\t|\tEBV\t|\t\t|\tacronym\t|\n"
    "10376\t|\tEpstein-Barr virus\t|\t\t|\tcommon name\t|\n"
)


@pytest.fixture
def cache_env(tmp_path, monkeypatch):
    """Redirect the taxdump cache under an isolated dir for this test."""
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    return tmp_path


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# parse_nodes_dmp
# ---------------------------------------------------------------------------


def test_parse_nodes_dmp_cache_hit_matches_cold(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    cold = F.parse_nodes_dmp(str(nodes))
    # Cold parse writes the cache.
    cache_path = F._taxdump_cache_path("nodes", str(nodes))
    assert cache_path.is_file()
    warm = F.parse_nodes_dmp(str(nodes))
    assert warm == cold
    assert warm[10376] == (3050299, "no rank")
    assert warm[3050299] == (3050294, "species")


def test_parse_nodes_dmp_reads_from_cache_not_source(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    F.parse_nodes_dmp(str(nodes))  # populate the cache
    # Poison the cache in place (same key, source unchanged). A genuine
    # cache hit must return the sentinel verbatim, proving the source is
    # not re-parsed.
    cache_path = F._taxdump_cache_path("nodes", str(nodes))
    sentinel = {999: (1, "sentinel")}
    F._cache_store(cache_path, sentinel)
    assert F.parse_nodes_dmp(str(nodes)) == sentinel


def test_parse_nodes_dmp_invalidates_on_source_change(cache_env):
    nodes = cache_env / "nodes.dmp"
    _write(nodes, NODES)
    first = F.parse_nodes_dmp(str(nodes))
    assert 200000 not in first
    # Append a node so both size and (after a beat) mtime differ; the
    # changed key misses the stale cache and re-parses.
    time.sleep(0.01)
    _write(nodes, NODES + "200000\t|\t10239\t|\tspecies\t|\n")
    second = F.parse_nodes_dmp(str(nodes))
    assert second[200000] == (10239, "species")


def test_parse_nodes_dmp_tolerates_corrupt_cache(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    cache_path = F._taxdump_cache_path("nodes", str(nodes))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(b"\x00 not a pickle \xff")
    # A garbage cache must not propagate; it falls back to a real parse.
    result = F.parse_nodes_dmp(str(nodes))
    assert result[10376] == (3050299, "no rank")
    # ... and the corrupt cache is replaced by a valid one.
    assert F._cache_load(cache_path) == result


def test_parse_nodes_dmp_missing_file_is_not_cached(cache_env):
    missing = cache_env / "absent.dmp"
    assert F.parse_nodes_dmp(str(missing)) == {}
    assert not F._taxdump_cache_path("nodes", str(missing)).is_file()


# ---------------------------------------------------------------------------
# _load_taxdump_for_species_walkup
# ---------------------------------------------------------------------------


def test_walkup_cache_hit_matches_cold(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    names = _write(cache_env / "names.dmp", NAMES)
    cold = F._load_taxdump_for_species_walkup(str(nodes), str(names))
    warm = F._load_taxdump_for_species_walkup(str(nodes), str(names))
    assert cold is not None and warm is not None
    node_info, sci_name, alias_name = warm
    assert (node_info, sci_name, alias_name) == cold
    assert sci_name[3050299] == "Lymphocryptovirus humangamma4"
    assert sci_name[10376] == "human gammaherpesvirus 4"
    assert "EBV" in alias_name[10376]
    assert "Epstein-Barr virus" in alias_name[10376]


def test_walkup_reads_from_cache_not_source(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    names = _write(cache_env / "names.dmp", NAMES)
    F._load_taxdump_for_species_walkup(str(nodes), str(names))
    cache_path = F._taxdump_cache_path("walkup", str(nodes), str(names))
    sentinel = ({1: (1, "root")}, {1: "sentinel sci"}, {1: ["sentinel alias"]})
    F._cache_store(cache_path, sentinel)
    assert F._load_taxdump_for_species_walkup(str(nodes), str(names)) == sentinel


def test_walkup_invalidates_on_names_change(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    names = cache_env / "names.dmp"
    _write(names, NAMES)
    first = F._load_taxdump_for_species_walkup(str(nodes), str(names))
    assert 55555 not in first[1]
    time.sleep(0.01)
    _write(names, NAMES + "55555\t|\tNovel virus\t|\t\t|\tscientific name\t|\n")
    second = F._load_taxdump_for_species_walkup(str(nodes), str(names))
    assert second[1][55555] == "Novel virus"


def test_walkup_missing_file_returns_none_uncached(cache_env):
    nodes = _write(cache_env / "nodes.dmp", NODES)
    absent = cache_env / "no_names.dmp"
    assert F._load_taxdump_for_species_walkup(str(nodes), str(absent)) is None
    assert not F._taxdump_cache_path("walkup", str(nodes), str(absent)).is_file()


# ---------------------------------------------------------------------------
# cache key + low-level helpers
# ---------------------------------------------------------------------------


def test_cache_path_separates_sources_kinds_and_tracks_size(cache_env):
    a = _write(cache_env / "a.dmp", NODES)
    b = _write(cache_env / "b.dmp", NODES + "9\t|\t1\t|\tno rank\t|\n")
    key_a = F._taxdump_cache_path("nodes", str(a))
    # Stable for the same kind + source.
    assert F._taxdump_cache_path("nodes", str(a)) == key_a
    # Distinct source -> distinct key.
    assert F._taxdump_cache_path("nodes", str(b)) != key_a
    # Distinct kind -> distinct key.
    assert F._taxdump_cache_path("walkup", str(a)) != key_a


def test_cache_path_changes_when_size_changes(cache_env):
    src = cache_env / "n.dmp"
    _write(src, NODES)
    before = F._taxdump_cache_path("nodes", str(src))
    _write(src, NODES + "9\t|\t1\t|\tno rank\t|\n")
    after = F._taxdump_cache_path("nodes", str(src))
    assert before != after


def test_cache_dir_is_user_namespaced(cache_env):
    # On a shared node, the cache subdir must be per-user so two users
    # cannot collide on (or be unable to overwrite) the same files.
    src = _write(cache_env / "n.dmp", NODES)
    path = F._taxdump_cache_path("nodes", str(src))
    assert path.parent.name == f"virushanter2_taxdump_{F._cache_user()}"


def test_cache_load_missing_returns_none(tmp_path):
    assert F._cache_load(tmp_path / "nope.pkl") is None


def test_cache_store_then_load_roundtrip(tmp_path):
    target = tmp_path / "nested" / "c.pkl"
    obj = {1: (2, "x"), 3: (4, "y")}
    F._cache_store(target, obj)  # creates the parent dir
    assert target.is_file()
    assert F._cache_load(target) == obj


def test_cache_store_is_silent_when_unwritable(tmp_path):
    # Parent path is a regular file, so the directory cannot be created;
    # the store must swallow the OSError rather than propagate it.
    blocker = tmp_path / "blocker"
    blocker.write_text("x")
    target = blocker / "child.pkl"
    F._cache_store(target, {1: 1})  # must not raise
    assert not target.exists()
