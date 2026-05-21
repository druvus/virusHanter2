"""Hermetic unit test for the rank filter + genus walk-up logic
that `bwa_align_to_kraken_hits` applies on every classifier hit.

The rule's `_record` closure cannot be imported directly from a
`.smk` file, so this test re-creates the same logic against
synthetic dictionaries and exercises the four expected branches:

1. Higher-rank taxids are dropped silently.
2. Species-level taxids present in the parquet are recorded.
3. Species-level taxids absent from the parquet but with a genus
   ancestor in the parquet are walked up; source tag carries the
   ``->genus`` suffix.
4. Species-level taxids absent from both the parquet and any
   genus ancestor in the parquet land in the unmapped sidecar.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_virus_parquet import find_genus_taxid  # noqa: E402


def _record_factory(
    parquet_tax_ids: set[int],
    nodes: dict[int, tuple[int, str]],
    rank_filter: set[str],
    coverage_genus_walkup: bool,
):
    sources_for_tid: dict[int, set[str]] = {}
    names_for_tid: dict[int, str] = {}
    unmapped_rows: list[tuple[int, str, str, str]] = []

    def _record(tid: int, name: str, source: str) -> None:
        rank = nodes.get(tid, (0, "unknown"))[1] if nodes else "unknown"
        if rank in rank_filter:
            return
        if tid in parquet_tax_ids:
            sources_for_tid.setdefault(tid, set()).add(source)
            if name and tid not in names_for_tid:
                names_for_tid[tid] = name
            return
        if coverage_genus_walkup and nodes:
            genus_tid = find_genus_taxid(tid, nodes)
            if genus_tid and genus_tid in parquet_tax_ids:
                sources_for_tid.setdefault(genus_tid, set()).add(
                    f"{source}->genus"
                )
                if name and genus_tid not in names_for_tid:
                    names_for_tid[genus_tid] = name
                return
        unmapped_rows.append((tid, name, source, "absent_from_parquet"))

    return _record, sources_for_tid, names_for_tid, unmapped_rows


# Synthetic NCBI-style mini-taxonomy:
#   1                                no rank (root)
#   10239                            superkingdom (Viruses)
#   687331                           genus  (Alphatorquevirus)
#   3048424                          species/strain (homin24)
#   3048433                          species/strain (homin9)
#   687329                           family (Anelloviridae)
#   10239 is also the parent of the family in this synthetic
#   chain; not realistic but irrelevant for the tests below.
_NODES: dict[int, tuple[int, str]] = {
    1: (1, "no rank"),
    10239: (1, "superkingdom"),
    687329: (10239, "family"),
    687331: (687329, "genus"),
    3048424: (687331, "species"),
    3048433: (687331, "species"),
}

_RANK_FILTER = {"family", "order", "class", "phylum", "kingdom", "realm"}


def test_rank_filtered_taxid_is_dropped_silently():
    record, sources, names, unmapped = _record_factory(
        parquet_tax_ids={687331},
        nodes=_NODES,
        rank_filter=_RANK_FILTER,
        coverage_genus_walkup=True,
    )
    # 687329 is a family — must be dropped before the parquet
    # lookup or the walk-up.
    record(687329, "Anelloviridae", "kraken")
    assert sources == {}
    assert unmapped == []


def test_species_in_parquet_is_recorded():
    record, sources, names, unmapped = _record_factory(
        parquet_tax_ids={3048424},
        nodes=_NODES,
        rank_filter=_RANK_FILTER,
        coverage_genus_walkup=True,
    )
    record(3048424, "Alphatorquevirus homin24", "kraken")
    assert sources == {3048424: {"kraken"}}
    assert names == {3048424: "Alphatorquevirus homin24"}
    assert unmapped == []


def test_species_walks_up_to_genus():
    """3048424 is a strain Kraken reported; absent from parquet;
    its genus 687331 *is* in parquet → walk-up substitutes it."""
    record, sources, names, unmapped = _record_factory(
        parquet_tax_ids={687331},
        nodes=_NODES,
        rank_filter=_RANK_FILTER,
        coverage_genus_walkup=True,
    )
    record(3048424, "Alphatorquevirus homin24", "kraken")
    assert sources == {687331: {"kraken->genus"}}
    assert names == {687331: "Alphatorquevirus homin24"}
    assert unmapped == []


def test_species_walkup_disabled_lands_in_unmapped():
    record, sources, names, unmapped = _record_factory(
        parquet_tax_ids={687331},
        nodes=_NODES,
        rank_filter=_RANK_FILTER,
        coverage_genus_walkup=False,
    )
    record(3048424, "Alphatorquevirus homin24", "kraken")
    assert sources == {}
    assert unmapped == [(3048424, "Alphatorquevirus homin24", "kraken", "absent_from_parquet")]


def test_species_no_genus_in_parquet_lands_in_unmapped():
    record, sources, names, unmapped = _record_factory(
        parquet_tax_ids=set(),  # nothing in parquet
        nodes=_NODES,
        rank_filter=_RANK_FILTER,
        coverage_genus_walkup=True,
    )
    record(3048424, "Alphatorquevirus homin24", "kraken")
    assert sources == {}
    assert unmapped == [(3048424, "Alphatorquevirus homin24", "kraken", "absent_from_parquet")]


def test_two_strains_collapse_to_one_genus_via_walkup():
    """Two strain-level Kraken hits with the same parent genus
    deduplicate at the genus level when the walk-up substitutes
    the same parquet reference for both. Both classifier names
    survive as alternative source tags."""
    record, sources, names, unmapped = _record_factory(
        parquet_tax_ids={687331},
        nodes=_NODES,
        rank_filter=_RANK_FILTER,
        coverage_genus_walkup=True,
    )
    record(3048424, "Alphatorquevirus homin24", "kraken")
    record(3048433, "Alphatorquevirus homin9", "kaiju")
    assert sources == {687331: {"kraken->genus", "kaiju->genus"}}
    # First non-empty name wins.
    assert names == {687331: "Alphatorquevirus homin24"}
