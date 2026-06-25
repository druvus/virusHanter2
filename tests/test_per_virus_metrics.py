"""Unit tests for the per-(sample, virus) join logic.

The metric extraction in `scripts/per_virus_metrics.py` is deterministic
once its inputs are framed as DataFrames + dicts, so the helpers can be
exercised without invoking Snakemake or any binary.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.per_virus_metrics import (  # noqa: E402
    attribute_contigs,
    build_per_virus_rows,
    kaiju_lookup,
    kraken_viral_top_n,
    mosdepth_summary_table,
    mosdepth_thresholds_table,
    parquet_first_token,
    parquet_refs_by_taxid,
    parse_fastp_total_reads,
    parse_flagstat,
    read_virus_parquet_taxids,
)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _kraken_df(rows):
    return pd.DataFrame(
        rows,
        columns=[
            "percent",
            "count_clades",
            "count",
            "tax_lvl",
            "taxonomy_id",
            "name",
            "domain",
        ],
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_parquet_first_token_strips_whitespace_tail():
    assert parquet_first_token("NC_000866.4 Enterobacteria phage T4") == "NC_000866.4"
    assert parquet_first_token("synthetic_virus") == "synthetic_virus"


def test_parse_flagstat_extracts_paired_and_mapped():
    text = (
        "400000 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "395102 + 0 primary\n"
        "0 + 0 secondary\n"
        "9799 + 0 supplementary\n"
        "0 + 0 duplicates\n"
        "395102 + 0 paired in sequencing\n"
        "26000 + 0 with itself and mate mapped\n"
    )
    total, mapped = parse_flagstat(text)
    assert total == 395102
    assert mapped == 26000


def test_parse_fastp_total_reads_reads_before_filtering():
    js = (
        '{"summary": {"before_filtering": {"total_reads": 1820,'
        ' "read1_mean_length": 150}}}'
    )
    assert parse_fastp_total_reads(js) == 1820


def test_kraken_viral_top_n_keeps_species_rank_rows_only():
    df = _kraken_df(
        [
            (10.0, 1000, 1000, "D", 10239, "Viruses", "Viruses"),
            (10.0, 1000, 1000, "F", 9999, "Some viral family", "Viruses"),
            (5.0, 500, 500, "S", 100001, "Phage X", "Viruses"),
            (60.0, 6000, 6000, "S", 100002, "Phage Y", "Viruses"),
            (45.0, 4500, 4500, "S1", 100004, "Phage Y subsp", "Viruses"),
            (40.0, 4000, 4000, "S2", 100005, "Phage Y strain", "Viruses"),
            (1.0, 100, 100, "S", 100003, "Phage Z", "Viruses"),
            (90.0, 9000, 9000, "D", 2, "Bacteria", "Bacteria"),
        ]
    )
    top = kraken_viral_top_n(df, top_n=3)
    # Sorted by percent desc, restricted to tax_lvl in {S, S1, S2} and
    # to the Viruses domain. Top 3: 60% S, 45% S1, 40% S2.
    # The 10% Viruses Domain (D) row and the 10% family (F) row are
    # filtered out by the species-rank guard; the bacterial row is
    # filtered out by the domain guard.
    assert list(top["taxonomy_id"]) == [100002, 100004, 100005]
    assert list(top["rank"]) == [0, 1, 2]


def test_kaiju_lookup_indexes_by_taxon_id_and_skips_na():
    df = pd.DataFrame(
        {
            "file": ["a", "a", "a"],
            "percent": [50.0, 0.0, 25.0],
            "reads": [100, 0, 50],
            "taxon_id": [100001, None, 100002],
            "taxon_name": ["Phage A", "unclassified", "Phage B"],
        }
    )
    out = kaiju_lookup(df)
    assert set(out.keys()) == {100001, 100002}
    assert out[100001]["reads"] == 100
    assert out[100002]["taxon_name"] == "Phage B"


def test_parquet_refs_by_taxid_groups_multiple_refs():
    df = pd.DataFrame(
        {
            "name": [
                "NC_001 strain 1",
                "NC_002 strain 2",
                "NC_003 strain only",
            ],
            "sequence": ["A", "C", "G"],
            "tax_id": [42, 42, 7],
        }
    )
    buckets = parquet_refs_by_taxid(df)
    assert set(buckets.keys()) == {42, 7}
    assert len(buckets[42]) == 2
    assert {ref["chrom"] for ref in buckets[42]} == {"NC_001", "NC_002"}
    assert buckets[42][0]["base_accession"] == "NC_001"


def test_mosdepth_summary_skips_total_rows():
    text = (
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        "synthetic_virus\t5000\t240000\t48.0\t0\t73\n"
        "total\t5000\t240000\t48.0\t0\t73\n"
        "total_region\t5000\t240000\t48.0\t0\t73\n"
    )
    out = mosdepth_summary_table(text)
    assert set(out.keys()) == {"synthetic_virus"}
    assert out["synthetic_virus"]["length"] == 5000


def test_mosdepth_thresholds_sums_per_chrom(tmp_path):
    bed_gz = tmp_path / "t.bed.gz"
    with gzip.open(bed_gz, "wt") as fh:
        fh.write("#chrom\tstart\tend\tregion\t1X\t5X\t10X\n")
        fh.write("syn_v\t0\t1000\tregion0\t1000\t800\t500\n")
        fh.write("syn_v\t1000\t2000\tregion1\t1000\t600\t100\n")
        fh.write("other_v\t0\t500\tregion0\t500\t250\t0\n")
    out = mosdepth_thresholds_table(bed_gz)
    assert out["syn_v"]["bases_ge_5x"] == 1400
    assert out["syn_v"]["bases_ge_1x"] == 2000
    assert out["other_v"]["bases_ge_10x"] == 0


def test_attribute_contigs_prefers_accession_match():
    blastn = pd.DataFrame(
        {
            "name": ["k21_0_pilon", "k21_1_pilon"],
            "match_name": ["irrelevant header text", "Murine reference seq"],
            "accession": ["NC_001.1", "NC_999.1"],
        }
    )
    parquet_acc_to_tax = {"NC_001": 42}
    # Use distinctive taxon names so the substring fallback is
    # unambiguous: only the second row's first token ("Murine") matches
    # taxid 7.
    taxid_to_name = {42: "Phage A", 7: "Murine virus"}
    contigs = attribute_contigs(blastn, parquet_acc_to_tax, taxid_to_name)
    # First row -> taxid 42 via accession (parquet hit), regardless of
    # the misleading match_name.
    assert contigs[42] == ["k21_0_pilon"]
    # Second row -> taxid 7 via substring of first token.
    assert contigs[7] == ["k21_1_pilon"]


def test_attribute_contigs_with_no_match_returns_empty():
    blastn = pd.DataFrame(
        {
            "name": ["k21_0_pilon"],
            "match_name": ["No such virus"],
            "accession": ["XX_000.0"],
        }
    )
    contigs = attribute_contigs(blastn, {}, {42: "Phage A"})
    assert contigs == {}


def test_build_per_virus_rows_aggregates_multi_reference_taxid(tmp_path):
    kraken = _kraken_df(
        [
            (10.0, 1000, 1000, "D", 10239, "Viruses", "Viruses"),
            (60.0, 600, 600, "S", 42, "Phage A", "Viruses"),
            (5.0, 50, 50, "S", 7, "Phage B", "Viruses"),
        ]
    )
    kaiju = pd.DataFrame(
        {
            "file": ["a", "a"],
            "percent": [80.0, 5.0],
            "reads": [600, 50],
            "taxon_id": [42, 7],
            "taxon_name": ["Phage A", "Phage B"],
        }
    )
    blastn = pd.DataFrame(
        {
            "match_name": ["Phage A strain 1", "Phage A strain 2", "Phage B refs"],
            "accession": ["NC_001.1", "NC_002.1", "NC_003.1"],
        }
    )
    parquet = pd.DataFrame(
        {
            "name": ["NC_001 strain 1", "NC_002 strain 2", "NC_003 only B"],
            "sequence": ["A", "C", "G"],
            "tax_id": [42, 42, 7],
        }
    )
    summary = {
        "NC_001": {"length": 1000, "bases": 4000, "mean": 4.0},
        "NC_002": {"length": 2000, "bases": 16000, "mean": 8.0},
        "NC_003": {"length": 500, "bases": 0, "mean": 0.0},
    }
    thresholds = {
        "NC_001": {"bases_ge_1x": 900, "bases_ge_5x": 200, "bases_ge_10x": 0},
        "NC_002": {"bases_ge_1x": 1800, "bases_ge_5x": 1500, "bases_ge_10x": 800},
        "NC_003": {"bases_ge_1x": 0, "bases_ge_5x": 0, "bases_ge_10x": 0},
    }
    df = build_per_virus_rows(
        run_name="251015_M00568_0723_000000000-DRRKK",
        sample_name="135_S1_L001_R",
        kraken_df=kraken,
        kaiju_df=kaiju,
        blastn_df=blastn,
        parquet_df=parquet,
        summary=summary,
        thresholds=thresholds,
        total_reads=10000,
        human_reads=2000,
        top_n=20,
    )
    # One row per kraken viral species (tax_lvl S/S1/S2). The
    # Domain "Viruses" row and any higher-rank ancestor rows are
    # filtered out so the per-virus CSV stays one-row-per-species.
    assert set(df["virus_taxid"]) == {42, 7}
    row_a = df.loc[df["virus_taxid"] == 42].iloc[0]
    row_b = df.loc[df["virus_taxid"] == 7].iloc[0]

    # Phage A spans NC_001 (length 1000) + NC_002 (length 2000) = 3000
    # bases of reference, 4000 + 16000 = 20000 bases aligned.
    assert row_a["mean_coverage"] == 20000 / 3000
    # bases_above_5x = 200 + 1500 = 1700; completeness = 1700/3000
    # rendered as a percent (0-100).
    assert row_a["bases_above_5x"] == 1700
    assert abs(row_a["Completeness (% >5X)"] - 100 * 1700 / 3000) < 1e-9
    # Both NC_001 and NC_002 BLASTN hits attribute to taxid 42.
    assert row_a["contigs"] == 2
    assert row_a["virus_name_kaiju"] == "Phage A"
    # Specific RPM uses count_clades=600 / 10000 * 1e6 = 60000.
    assert row_a["specific_virus_rpm"] == 60000.0
    # All-virus RPM uses the Domain row's count_clades (1000) which
    # already includes every clade beneath it. 1000 / 10000 * 1e6 = 100000.
    assert row_a["all_virus_rpm"] == 100000.0
    assert row_b["all_virus_rpm"] == 100000.0
    # Other reads = (10000 - 2000) - 1000 = 7000.
    assert row_a["other_reads"] == 7000
    # Per-row constants.
    assert row_a["total_reads"] == 10000
    assert row_a["human_reads"] == 2000
    assert row_a["human_reads_percent"] == 20.0
    assert row_a["non_human_reads"] == 8000
    assert row_a["non_human_reads_percent"] == 80.0
    assert row_a["note"] == ""
    assert row_a["run_name"] == "251015_M00568_0723_000000000-DRRKK"
    assert row_a["date"] == "251015"


def test_build_per_virus_rows_viral_only_kraken_uses_r1_anchor():
    """Mirror the smaller k2_viral_* DB shape where Viruses sits at
    tax_lvl == 'R1'. The Domain row's count_clades must still feed
    all_viral_reads and the resulting all_virus_rpm.
    """
    kraken = _kraken_df(
        [
            (3.0, 30, 30, "U", 0, "unclassified", "unclassified"),
            (97.0, 970, 0, "R", 1, "root", "root"),
            (97.0, 970, 0, "R1", 10239, "Viruses", "Viruses"),
            (60.0, 600, 600, "S", 42, "Phage A", "Viruses"),
        ]
    )
    df = build_per_virus_rows(
        run_name="251015_test",
        sample_name="s_R",
        kraken_df=kraken,
        kaiju_df=pd.DataFrame(),
        blastn_df=pd.DataFrame(),
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={},
        thresholds={},
        total_reads=1000,
        human_reads=30,
        top_n=10,
    )
    # The R1 anchor is used to derive all_viral_reads (970), but the
    # R1 row itself does not appear in the per-virus DataFrame
    # because the species-rank filter rejects R1. Only Phage A (S)
    # survives the filter; its all_virus_rpm must still reflect the
    # 970 reads from the R1 row's count_clades (not 0 as it would
    # be with a D-only anchor).
    row = df.loc[df["virus_name_kraken"] == "Phage A"].iloc[0]
    assert row["all_virus_rpm"] == 970 * 1_000_000.0 / 1000
    assert row["other_reads"] == (1000 - 30) - 970


def test_build_per_virus_rows_flags_dummy_contig_in_note():
    """When the BLASTN merged CSV's only contig is the DUMMY_CONTIG
    Pilon-polished placeholder, every per-virus row should carry a
    note flagging the silent MEGAHIT failure so reviewers do not read
    the report as a clean negative result.
    """
    kraken = _kraken_df(
        [
            (10.0, 100, 100, "D", 10239, "Viruses", "Viruses"),
            (5.0, 50, 50, "S", 42, "Phage A", "Viruses"),
        ]
    )
    # bblastn merged CSV: only the dummy contig made it through.
    blastn = pd.DataFrame(
        {
            "name": ["DUMMY_CONTIG_pilon"],
            "match_name": [""],
            "accession": [""],
        }
    )
    df = build_per_virus_rows(
        run_name="x",
        sample_name="s",
        kraken_df=kraken,
        kaiju_df=pd.DataFrame(),
        blastn_df=blastn,
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={},
        thresholds={},
        total_reads=1000,
        human_reads=10,
        top_n=10,
    )
    assert not df.empty
    assert (df["note"] == "MEGAHIT assembly failed; dummy contig only").all()


def test_dummy_sentinel_round_trips_to_assembly_failure_note(tmp_path):
    """End-to-end producer -> CSV -> consumer check.

    `dummy_contig_sentinel` (the merge_checkv_blastn producer) must emit
    a row that, once written to the merged CSV and read back the way
    per_virus_metrics.main does, drives `build_per_virus_rows` to flag
    the assembly failure. This guards the NaN-filled sentinel schema
    against breaking attribution downstream.
    """
    from scripts.functions import dummy_contig_sentinel

    merged_cols = ["name", "assembler", "match_name", "accession", "tax_id"]
    merged = pd.DataFrame(columns=merged_cols)  # empty inner join
    checkv = pd.DataFrame({"name": ["DUMMY_CONTIG_pilon"], "viral_genes": [0]})

    merged = dummy_contig_sentinel(merged, checkv, "MEGAHIT")
    csv_path = tmp_path / "merged.csv"
    merged.to_csv(csv_path, index=False)

    # Read back exactly as per_virus_metrics.main does.
    blastn_df = pd.read_csv(csv_path)

    kraken = _kraken_df(
        [
            (10.0, 100, 100, "D", 10239, "Viruses", "Viruses"),
            (5.0, 50, 50, "S", 42, "Phage A", "Viruses"),
        ]
    )
    df = build_per_virus_rows(
        run_name="x",
        sample_name="s",
        kraken_df=kraken,
        kaiju_df=pd.DataFrame(),
        blastn_df=blastn_df,
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={},
        thresholds={},
        total_reads=1000,
        human_reads=10,
        top_n=10,
    )
    assert not df.empty
    assert (df["note"] == "MEGAHIT assembly failed; dummy contig only").all()
    # The dummy must not be attributed as a real contig to any virus.
    assert (df["contigs"] == 0).all()


def test_build_per_virus_rows_dedupes_by_virus_name_keeping_max_contigs():
    """Kraken's S / S1 / S2 sub-rank chain often produces multiple
    rows that all canonicalise to the same virus name. The output
    should keep one row per virus_name_kraken - the one carrying
    the most attributed contigs (typically the S1 leaf where the
    parquet reference accession lives)."""
    kraken = _kraken_df(
        [
            # ICTV root (S) - no parquet reference, no contigs.
            (60.0, 600, 0, "S", 3050293, "Phage A", "Viruses"),
            # Legacy NCBI name (S1) - parquet keys here; contigs land here.
            (60.0, 600, 600, "S1", 42, "Phage A", "Viruses"),
            # Strain (S2) - no parquet reference.
            (5.0, 50, 50, "S2", 99999, "Phage A", "Viruses"),
            # Distinct virus, single row.
            (10.0, 100, 100, "S", 7, "Phage B", "Viruses"),
        ]
    )
    blastn = pd.DataFrame(
        {
            "name": ["k57_0_pilon", "k57_1_pilon", "k57_2_pilon"],
            "match_name": ["Phage A", "Phage A", "Phage A"],
            "accession": ["NC_001.1", "NC_001.1", "NC_001.1"],
        }
    )
    df = build_per_virus_rows(
        run_name="x",
        sample_name="s",
        kraken_df=kraken,
        kaiju_df=pd.DataFrame(),
        blastn_df=blastn,
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={},
        thresholds={},
        total_reads=1000,
        human_reads=0,
        top_n=10,
    )
    # Two viruses survive: Phage A and Phage B. Phage A's three
    # taxonomic rows collapse to the S1 row with the 3 contigs.
    assert len(df) == 2
    phage_a = df.loc[df["virus_name_kraken"] == "Phage A"].iloc[0]
    assert phage_a["virus_taxid"] == 42
    assert phage_a["contigs"] == 3


def _build_minimal(
    kraken_rows,
    *,
    kaiju_df=None,
    blastn_df=None,
    parquet_df=None,
    summary=None,
    thresholds=None,
    total_reads=1000,
    human_reads=0,
):
    """Boilerplate-free wrapper around ``build_per_virus_rows`` for the
    note-rule tests. Each rule only needs a tiny subset of the full
    input set, so the defaults are empty frames / empty dicts."""
    return build_per_virus_rows(
        run_name="x",
        sample_name="s",
        kraken_df=_kraken_df(kraken_rows),
        kaiju_df=kaiju_df if kaiju_df is not None else pd.DataFrame(),
        blastn_df=blastn_df if blastn_df is not None else pd.DataFrame(),
        parquet_df=parquet_df if parquet_df is not None else pd.DataFrame(),
        summary=summary or {},
        thresholds=thresholds or {},
        total_reads=total_reads,
        human_reads=human_reads,
        top_n=10,
    )


@pytest.mark.parametrize(
    "kraken_rows,expected_note",
    [
        pytest.param(
            [(5.0, 7, 7, "S", 42, "Phage A", "Viruses")],
            "low Kraken support (<10 reads)",
            id="low-kraken-support",
        ),
        pytest.param(
            [(50.0, 500, 500, "S", 42, "Phage A", "Viruses")],
            "no contigs attributed despite >=100 Kraken reads",
            id="no-contigs-despite-reads",
        ),
    ],
)
def test_per_row_note_simple_rules(kraken_rows, expected_note):
    """Note rules that need only a Kraken row to fire."""
    df = _build_minimal(kraken_rows)
    assert df.iloc[0]["note"] == expected_note


def test_per_row_note_flags_kraken_kaiju_disagreement():
    """Both classifiers have a species call but the names differ
    after normalisation - flag the disagreement."""
    df = _build_minimal(
        [(50.0, 500, 500, "S", 42, "Simplexvirus paninealpha3", "Viruses")],
        kaiju_df=pd.DataFrame(
            {
                "file": ["a"],
                "percent": [50.0],
                "reads": [500],
                "taxon_id": [42],
                "taxon_name": ["Simplexvirus humanalpha2"],
            }
        ),
    )
    assert "Kraken/Kaiju species disagree" in df.iloc[0]["note"]


def test_per_row_note_flags_low_mean_coverage():
    """Mean coverage strictly between 0 and 1x should fire the
    'mean coverage <1x' rule. Reference length zero is silent (no
    parquet reference means the rule cannot meaningfully fire)."""
    df = _build_minimal(
        [(50.0, 500, 500, "S", 42, "Phage A", "Viruses")],
        blastn_df=pd.DataFrame(
            {
                "name": ["k57_0_pilon"] * 200,
                "match_name": ["Phage A"] * 200,
                "accession": ["NC_001.1"] * 200,
            }
        ),
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={"NC_001": {"length": 1000, "bases": 500, "mean": 0.5}},
        thresholds={"NC_001": {"bases_ge_1x": 100, "bases_ge_5x": 0, "bases_ge_10x": 0}},
    )
    assert "mean coverage <1x" in df.iloc[0]["note"]


def test_build_per_virus_rows_real_contigs_leave_note_empty():
    """Sanity check the opposite branch — a real contig keeps note=''."""
    kraken = _kraken_df(
        [
            (10.0, 100, 100, "D", 10239, "Viruses", "Viruses"),
            (5.0, 50, 50, "S", 42, "Phage A", "Viruses"),
        ]
    )
    blastn = pd.DataFrame(
        {
            "name": ["k57_0_pilon"],
            "match_name": ["Phage A"],
            "accession": ["NC_001.1"],
        }
    )
    df = build_per_virus_rows(
        run_name="x",
        sample_name="s",
        kraken_df=kraken,
        kaiju_df=pd.DataFrame(),
        blastn_df=blastn,
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={},
        thresholds={},
        total_reads=1000,
        human_reads=10,
        top_n=10,
    )
    assert (df["note"] == "").all()


def test_build_per_virus_rows_zero_total_reads_does_not_divide_by_zero():
    kraken = _kraken_df(
        [
            (10.0, 0, 0, "D", 10239, "Viruses", "Viruses"),
            (5.0, 0, 0, "S", 42, "Phage A", "Viruses"),
        ]
    )
    df = build_per_virus_rows(
        run_name="x",
        sample_name="s",
        kraken_df=kraken,
        kaiju_df=pd.DataFrame(),
        blastn_df=pd.DataFrame(),
        parquet_df=pd.DataFrame(
            {"name": ["NC_001 a"], "sequence": ["A"], "tax_id": [42]}
        ),
        summary={},
        thresholds={},
        total_reads=0,
        human_reads=0,
        top_n=10,
    )
    row = df.iloc[0]
    assert row["specific_virus_rpm"] == 0.0
    assert row["all_virus_rpm"] == 0.0
    assert row["human_reads_percent"] == 0.0


# ---------------------------------------------------------------------------
# read_virus_parquet_taxids: column projection + fallback
# ---------------------------------------------------------------------------


def _write_parquet(path: Path, frame: pd.DataFrame) -> Path:
    frame.to_parquet(path)
    return path


def test_read_virus_parquet_taxids_projects_only_needed_columns(tmp_path):
    # A full-schema parquet; the reader must skip the heavy `sequence`
    # column and return name + tax_id only.
    p = _write_parquet(
        tmp_path / "v.parquet",
        pd.DataFrame(
            {
                "name": ["NC_001 alpha", "NC_002 beta"],
                "sequence": ["ACGT" * 1000, "TTTT" * 1000],
                "tax_id": [100001, 100002],
                "rank": ["species", "species"],
                "genus_taxid": [1, 2],
            }
        ),
    )
    out = read_virus_parquet_taxids(p)
    assert list(out.columns) == ["name", "tax_id"]
    assert "sequence" not in out.columns
    assert out["tax_id"].tolist() == [100001, 100002]


def test_read_virus_parquet_taxids_feeds_parquet_refs_by_taxid(tmp_path):
    # The projected frame must still drive the downstream bucketer that
    # only reads name + tax_id, proving the projection is sufficient.
    p = _write_parquet(
        tmp_path / "v.parquet",
        pd.DataFrame(
            {
                "name": ["NC_001.1 alpha", "NC_009.2 beta"],
                "sequence": ["ACGT", "TTTT"],
                "tax_id": [100001, 100002],
            }
        ),
    )
    buckets = parquet_refs_by_taxid(read_virus_parquet_taxids(p))
    assert buckets[100001][0]["chrom"] == "NC_001.1"
    assert buckets[100001][0]["base_accession"] == "NC_001"


def test_read_virus_parquet_taxids_falls_back_when_column_absent(tmp_path):
    # An older parquet build without `tax_id`: pyarrow raises
    # ArrowInvalid (a ValueError subclass) for the missing projected
    # field, and the reader must fall back to a full read instead of
    # propagating the error.
    p = _write_parquet(
        tmp_path / "old.parquet",
        pd.DataFrame({"name": ["NC_001 alpha"], "sequence": ["ACGT"]}),
    )
    out = read_virus_parquet_taxids(p)
    assert "name" in out.columns  # did not raise; fell back to full read
    assert "sequence" in out.columns
