"""Extract per-(sample, virus) metrics for the Twist VRP collaborator.

Joins the existing pipeline outputs (Kraken2 report, Kaiju table,
BLASTN merged CSV, mosdepth summary + thresholds, fastp JSON, host
flagstat) and the workflow-level viral parquet into a flat CSV with
one row per detected virus (Kraken taxid) per sample.

The schema is described in `docs/PER_VIRUS_OUTPUT.md`.

Usage:

    python scripts/per_virus_metrics.py \\
        --sample-name 135_S1_L001_R \\
        --run-name 251015_M00568_0723_000000000-DRRKK \\
        --kraken-csv      .../KRAKEN/135.kraken.csv \\
        --kaiju-tsv       .../KAIJU/135.kaiju.table.tsv \\
        --blastn-csv      .../CHECKV/135.merged.csv \\
        --mosdepth-summary    .../MOSDEPTH/135.mosdepth.summary.txt \\
        --mosdepth-thresholds .../MOSDEPTH/135.thresholds.bed.gz \\
        --fastp-json      .../FASTP/135.fastp.json \\
        --flagstat        .../logs/human_contamination_flagstat.txt \\
        --virus-parquet   .../INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet \\
        --top-n 20 \\
        --out             .../135.per_virus.csv
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Pure helpers — unit-tested in tests/test_per_virus_metrics.py.
# ---------------------------------------------------------------------------


_FLAGSTAT_TOTAL_RE = re.compile(r"^(\d+)\s+\+\s+\d+\s+paired in sequencing", re.M)
_FLAGSTAT_MAPPED_RE = re.compile(
    r"^(\d+)\s+\+\s+\d+\s+with itself and mate mapped", re.M
)


def parse_flagstat(text: str) -> tuple[int, int]:
    """Return (paired_in_sequencing, with_itself_and_mate_mapped) from a
    `samtools flagstat` report.
    """
    total_match = _FLAGSTAT_TOTAL_RE.search(text)
    mapped_match = _FLAGSTAT_MAPPED_RE.search(text)
    total = int(total_match.group(1)) if total_match else 0
    mapped = int(mapped_match.group(1)) if mapped_match else 0
    return total, mapped


def parse_fastp_total_reads(fastp_json_text: str) -> int:
    """Return `summary.before_filtering.total_reads` from a fastp JSON."""
    doc = json.loads(fastp_json_text)
    return int(doc.get("summary", {}).get("before_filtering", {}).get("total_reads", 0))


def kraken_viral_top_n(kraken_df: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    """Replicate `bwa_align_to_kraken_hits`'s selection: take the rows in
    the Viruses domain, sort by percent descending, keep the top N.

    Returns the same columns as the input, plus a stable `rank` 0..N-1.
    """
    viral = (
        kraken_df.loc[kraken_df["domain"] == "Viruses"]
        .sort_values("percent", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    viral["rank"] = viral.index
    return viral


def parquet_first_token(name: str) -> str:
    """First whitespace-delimited token of a FASTA header, matching the
    chrom name BWA / mosdepth produce.
    """
    return str(name).split()[0]


def parquet_base_accession(name: str) -> str:
    """First token without the trailing `.VERSION` suffix.
    `NC_000866.4 ...` -> `NC_000866`.
    """
    return parquet_first_token(name).split(".")[0]


def kaiju_lookup(kaiju_df: pd.DataFrame) -> dict[int, dict[str, object]]:
    """Index a Kaiju TSV by integer `taxon_id`. NA taxon_ids are dropped
    (those are the "unclassified" / "cannot be assigned" rows).
    """
    if kaiju_df.empty or "taxon_id" not in kaiju_df.columns:
        return {}
    tid = pd.to_numeric(kaiju_df["taxon_id"], errors="coerce")
    df = kaiju_df.assign(tid_int=tid).dropna(subset=["tid_int"])
    out: dict[int, dict[str, object]] = {}
    for row in df.itertuples():
        out[int(row.tid_int)] = {
            "taxon_name": getattr(row, "taxon_name", ""),
            "reads": int(getattr(row, "reads", 0) or 0),
            "percent": float(getattr(row, "percent", 0.0) or 0.0),
        }
    return out


def parquet_refs_by_taxid(parquet_df: pd.DataFrame) -> dict[int, list[dict]]:
    """Bucket parquet rows by tax_id. Each bucket holds one dict per
    reference with the first-token chrom id and the base accession (so
    BLAST hits without an explicit version still match).
    """
    out: dict[int, list[dict]] = {}
    for row in parquet_df.itertuples():
        try:
            tid = int(row.tax_id)
        except (ValueError, TypeError):
            continue
        if tid == 0:
            continue
        chrom = parquet_first_token(row.name)
        base = parquet_base_accession(row.name)
        out.setdefault(tid, []).append({"chrom": chrom, "base_accession": base})
    return out


def mosdepth_summary_table(summary_text: str) -> dict[str, dict[str, float]]:
    """Parse `mosdepth.summary.txt` (TSV with a header) into
    `{chrom: {length, bases, mean}}`. The trailing `total` and
    `total_region` rows are skipped.
    """
    out: dict[str, dict[str, float]] = {}
    for i, line in enumerate(summary_text.splitlines()):
        if i == 0 or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        chrom = parts[0]
        if chrom in ("total", "total_region"):
            continue
        try:
            length = int(parts[1])
            bases = int(parts[2])
            mean = float(parts[3])
        except ValueError:
            continue
        out[chrom] = {"length": length, "bases": bases, "mean": mean}
    return out


def mosdepth_thresholds_table(thresholds_bed_gz: Path) -> dict[str, dict[str, int]]:
    """Sum the per-region 1X/5X/10X columns of a mosdepth thresholds
    BED into `{chrom: {bases_>=1, bases_>=5, bases_>=10}}`.

    The file ships with a `#chrom  start  end  region  1X  5X  10X`
    header; subsequent rows are TSV.
    """
    out: dict[str, dict[str, int]] = {}
    header: list[str] | None = None
    with gzip.open(thresholds_bed_gz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [p.lstrip("#") for p in parts]
                continue
            if header is None or len(parts) < len(header):
                continue
            chrom = parts[0]
            row_counts = out.setdefault(
                chrom, {"bases_ge_1x": 0, "bases_ge_5x": 0, "bases_ge_10x": 0}
            )
            for col_name, raw in zip(header, parts):
                if col_name == "1X":
                    row_counts["bases_ge_1x"] += int(raw or 0)
                elif col_name == "5X":
                    row_counts["bases_ge_5x"] += int(raw or 0)
                elif col_name == "10X":
                    row_counts["bases_ge_10x"] += int(raw or 0)
    return out


def attribute_contigs(
    blastn_df: pd.DataFrame,
    parquet_acc_to_tax: dict[str, int],
    taxid_to_name: dict[int, str],
) -> dict[int, list[str]]:
    """Group BLASTN-merged rows by the Kraken taxid each is attributable to.

    Two-tier match: by accession (preferred), then by case-insensitive
    substring of the BLASTN `match_name`'s first token in the Kraken
    name. Returns ``{taxid: [contig_name, ...]}``; per-taxid contig
    counts are then ``len(values)``.
    """
    by_taxid: dict[int, list[str]] = {}
    if blastn_df.empty:
        return by_taxid
    for row in blastn_df.itertuples():
        tid: int | None = None
        # Tier 1: accession lookup.
        accession = getattr(row, "accession", None)
        if accession is not None and not pd.isna(accession):
            base = str(accession).split(".")[0]
            tid = parquet_acc_to_tax.get(base) or parquet_acc_to_tax.get(
                str(accession)
            )
        # Tier 2: substring of match_name in any kraken name.
        if tid is None:
            match_name = getattr(row, "match_name", None)
            if match_name is not None and not pd.isna(match_name):
                tokens = str(match_name).split()
                hit_token = tokens[0].lower() if tokens else ""
                for candidate_tid, kraken_name in taxid_to_name.items():
                    if hit_token and hit_token in str(kraken_name).lower():
                        tid = candidate_tid
                        break
        if tid is not None:
            contig = getattr(row, "name", None)
            by_taxid.setdefault(tid, []).append(
                "" if contig is None or pd.isna(contig) else str(contig)
            )
    return by_taxid


def parse_genomad_summary(path: Path) -> dict[str, float]:
    """Parse a geNomad ``<sample>_virus_summary.tsv`` into a per-contig
    map of ``{seq_name_first_token: virus_score}``.

    geNomad scores contigs and reports a virus_score in [0, 1]; rows that
    do not pass geNomad's own viral threshold are absent from this file,
    so any contig in the dict is a geNomad-called virus. Missing or
    unparseable scores are dropped silently.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}
    df = pd.read_csv(path, sep="\t")
    if "seq_name" not in df.columns or "virus_score" not in df.columns:
        return {}
    out: dict[str, float] = {}
    for row in df.itertuples():
        try:
            seq = str(row.seq_name).split()[0]
            score = float(row.virus_score)
        except (AttributeError, TypeError, ValueError):
            continue
        out[seq] = score
    return out


def build_per_virus_rows(
    *,
    run_name: str,
    sample_name: str,
    kraken_df: pd.DataFrame,
    kaiju_df: pd.DataFrame,
    blastn_df: pd.DataFrame,
    parquet_df: pd.DataFrame,
    summary: dict[str, dict[str, float]],
    thresholds: dict[str, dict[str, int]],
    total_reads: int,
    human_reads: int,
    top_n: int = 20,
    genomad_summary: dict[str, float] | None = None,
    assemblers: list[str] | None = None,
) -> pd.DataFrame:
    """Top-level join. Returns a DataFrame with the schema documented in
    `docs/PER_VIRUS_OUTPUT.md`.
    """
    viral = kraken_viral_top_n(kraken_df, top_n=top_n)

    # All-virus reads = the Domain "Viruses" row's `count_clades`,
    # which already accounts for every species clade beneath it.
    # Summing the whole `domain == "Viruses"` subtree would
    # double-count (Domain + all descendants). Accept either tax_lvl
    # "D" (pluspf) or "R1" (viral-only DBs such as k2_viral_*).
    # Matches the `kraken_virus_percent` parity fix in
    # `aggregate_run_information.py`.
    domain_rows = kraken_df.loc[
        (kraken_df["tax_lvl"].isin(["D", "R1"])) & (kraken_df["name"] == "Viruses"),
        "count_clades",
    ]
    all_viral_reads = int(domain_rows.iloc[0]) if not domain_rows.empty else 0
    non_human_reads = total_reads - human_reads
    human_pct = 100.0 * human_reads / total_reads if total_reads > 0 else 0.0
    non_human_pct = 100.0 - human_pct if total_reads > 0 else 0.0
    other_reads = non_human_reads - all_viral_reads
    all_virus_rpm = (
        all_viral_reads * 1_000_000.0 / total_reads if total_reads > 0 else 0.0
    )

    parquet_buckets = parquet_refs_by_taxid(parquet_df)
    parquet_acc_to_tax: dict[str, int] = {}
    for tid, refs in parquet_buckets.items():
        for ref in refs:
            parquet_acc_to_tax[ref["base_accession"]] = tid

    kaiju_by_tid = kaiju_lookup(kaiju_df)
    taxid_to_name = {int(r.taxonomy_id): r.name for r in viral.itertuples()}
    contigs_by_taxid = attribute_contigs(blastn_df, parquet_acc_to_tax, taxid_to_name)

    # Map contig name -> assembler so we can split per-taxid contig
    # counts by assembler without re-walking the BLAST CSV. When
    # `assembler` is absent (single-assembler run, no wrangle_pilon
    # touch yet, or a legacy CSV) every contig is attributed to a
    # blank assembler and the per-assembler columns are all zero.
    contig_to_assembler: dict[str, str] = {}
    if not blastn_df.empty and "assembler" in blastn_df.columns:
        for r in blastn_df.itertuples():
            name = getattr(r, "name", None)
            asm = getattr(r, "assembler", "") or ""
            if name is not None and not pd.isna(name):
                contig_to_assembler[str(name)] = str(asm)

    # Detect the MEGAHIT-failure-fallback case: the upstream rule
    # writes a single DUMMY_CONTIG sequence when MEGAHIT crashes (most
    # commonly on Apple Silicon at small k). Pilon then polishes that
    # into ``DUMMY_CONTIG_pilon``. Without surfacing it, downstream
    # per-virus rows report "contigs: 0" with otherwise normal reads
    # counts and look indistinguishable from a low-coverage real run.
    # Flag every per-virus row of such samples so reviewers can spot
    # silent assembly failures.
    sample_note = ""
    if not blastn_df.empty and "name" in blastn_df.columns:
        contig_names = blastn_df["name"].dropna().astype(str)
        if (
            len(contig_names) > 0
            and contig_names.str.startswith("DUMMY_CONTIG").all()
        ):
            sample_note = "MEGAHIT assembly failed; dummy contig only"

    date_part = run_name.split("_")[0] if run_name else ""

    rows: list[dict] = []
    for v in viral.itertuples():
        taxid = int(v.taxonomy_id)
        virus_reads_kraken = int(v.count_clades)
        # Aggregate mosdepth across all references for this taxid.
        total_length = 0
        total_bases = 0
        total_ge5x = 0
        for ref in parquet_buckets.get(taxid, []):
            chrom = ref["chrom"]
            s = summary.get(chrom)
            if s is None:
                continue
            total_length += int(s["length"])
            total_bases += int(s["bases"])
            t = thresholds.get(chrom, {})
            total_ge5x += int(t.get("bases_ge_5x", 0))
        mean_cov = total_bases / total_length if total_length > 0 else 0.0
        completeness_5x = total_ge5x / total_length if total_length > 0 else 0.0

        specific_rpm = (
            virus_reads_kraken * 1_000_000.0 / total_reads if total_reads > 0 else 0.0
        )

        kaiju_row = kaiju_by_tid.get(taxid, {})
        attributed_contigs = contigs_by_taxid.get(taxid, [])

        # geNomad columns are additive trailing columns: present only
        # when a geNomad summary was supplied. Counts the attributed
        # contigs that geNomad called viral; ``max_score`` is the
        # highest virus_score among them (NaN when none qualify).
        row = {
            "run_name": run_name,
            "sample_name": sample_name,
            "date": date_part,
            "virus_name_kraken": v.name,
            "virus_taxid": taxid,
            "virus_name_kaiju": kaiju_row.get("taxon_name", ""),
            "contigs": len(attributed_contigs),
            "virus_reads_kraken2": virus_reads_kraken,
            "other_reads": other_reads,
            "total_reads": total_reads,
            "human_reads": human_reads,
            "human_reads_percent": human_pct,
            "non_human_reads": non_human_reads,
            "non_human_reads_percent": non_human_pct,
            "note": sample_note,
            "specific_virus_rpm": specific_rpm,
            "all_virus_rpm": all_virus_rpm,
            "completeness_5x": completeness_5x,
            "bases_above_5x": total_ge5x,
            "mean_coverage": mean_cov,
        }
        # Per-assembler contig counts. Trailing columns; absent when
        # the caller did not declare an `assemblers` list (legacy
        # single-assembler call sites stay byte-identical).
        if assemblers:
            for asm in assemblers:
                col = f"{asm.lower()}_contigs"
                row[col] = sum(
                    1
                    for c in attributed_contigs
                    if contig_to_assembler.get(c, "") == asm
                )
        if genomad_summary is not None:
            matched = [
                genomad_summary[c]
                for c in attributed_contigs
                if c in genomad_summary
            ]
            row["genomad_viral_contigs"] = len(matched)
            row["genomad_max_virus_score"] = max(matched) if matched else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CLI entry point.
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--sample-name", required=True)
    p.add_argument("--run-name", required=True)
    p.add_argument("--kraken-csv", required=True, type=Path)
    p.add_argument("--kaiju-tsv", required=True, type=Path)
    p.add_argument(
        "--blastn-csv",
        required=True,
        type=Path,
        action="append",
        help=(
            "Per-assembler BLASTN merged CSV. Repeat the flag once per "
            "assembler; the CSVs are concatenated and the `assembler` "
            "column carries through to the per-virus rows."
        ),
    )
    p.add_argument("--mosdepth-summary", required=True, type=Path)
    p.add_argument("--mosdepth-thresholds", required=True, type=Path)
    p.add_argument("--fastp-json", required=True, type=Path)
    p.add_argument("--flagstat", required=True, type=Path)
    p.add_argument("--virus-parquet", required=True, type=Path)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument(
        "--genomad-summary",
        type=Path,
        default=None,
        action="append",
        help=(
            "Optional geNomad <sample>_virus_summary.tsv. Repeat once "
            "per assembler when geNomad is enabled. When supplied, "
            "two trailing columns (genomad_viral_contigs and "
            "genomad_max_virus_score) are appended to the per-virus CSV."
        ),
    )
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    kraken_df = pd.read_csv(args.kraken_csv)
    kaiju_df = (
        pd.read_csv(args.kaiju_tsv, sep="\t") if args.kaiju_tsv.exists() else pd.DataFrame()
    )

    # Concatenate every per-assembler BLAST/CheckV merged CSV. Each
    # has an `assembler` column written by `wrangle_pilon`. Skip
    # missing or empty files quietly so a failed assembler does not
    # tear the per-virus aggregation down.
    blastn_frames: list[pd.DataFrame] = []
    declared_assemblers: list[str] = []
    for csv in args.blastn_csv:
        if not csv.exists() or csv.stat().st_size == 0:
            continue
        df_part = pd.read_csv(csv)
        if "assembler" not in df_part.columns:
            df_part = df_part.assign(assembler="")
        blastn_frames.append(df_part)
        for asm in df_part["assembler"].dropna().unique().tolist():
            if asm and asm not in declared_assemblers:
                declared_assemblers.append(str(asm))
    blastn_df = (
        pd.concat(blastn_frames, ignore_index=True)
        if blastn_frames
        else pd.DataFrame()
    )

    parquet_df = pd.read_parquet(args.virus_parquet)

    summary = mosdepth_summary_table(args.mosdepth_summary.read_text())
    thresholds = mosdepth_thresholds_table(args.mosdepth_thresholds)

    total_reads = parse_fastp_total_reads(args.fastp_json.read_text())
    _, human_reads = parse_flagstat(args.flagstat.read_text())

    # Merge per-assembler geNomad summaries into a single contig ->
    # score dict; identical contig names across assemblers are
    # vanishingly unlikely (each carries its own k-mer / NODE name)
    # so a flat dict is sufficient.
    genomad_summary: dict[str, float] | None
    if args.genomad_summary:
        merged: dict[str, float] = {}
        for p in args.genomad_summary:
            merged.update(parse_genomad_summary(p))
        genomad_summary = merged
    else:
        genomad_summary = None

    df = build_per_virus_rows(
        run_name=args.run_name,
        sample_name=args.sample_name,
        kraken_df=kraken_df,
        kaiju_df=kaiju_df,
        blastn_df=blastn_df,
        parquet_df=parquet_df,
        summary=summary,
        thresholds=thresholds,
        total_reads=total_reads,
        human_reads=human_reads,
        top_n=args.top_n,
        genomad_summary=genomad_summary,
        assemblers=declared_assemblers if declared_assemblers else None,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
