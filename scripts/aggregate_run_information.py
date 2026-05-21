"""Aggregate per-sample metrics into a run-level CSV.

Invoked by the `aggregate_run_information` Snakemake rule via the `script:`
directive, so it runs inside the reporthanter conda env. The Snakemake
object exposes `snakemake.input`, `snakemake.output`, `snakemake.params`,
and `snakemake.config`.

Note: a `from __future__ import annotations` would not work here because
Snakemake prepends its own boilerplate when materialising script: rules.
"""
import json
from pathlib import Path

import pandas as pd

from reporthanter import FlagstatProcessor

# Pull in the pipeline-side helper for hex-encoding the HTML blob.
from scripts.functions import read_file_as_blob


def aggregate_sample_info(sample_folder: Path) -> pd.DataFrame:
    sample_folder = sample_folder.resolve()
    sample_name = sample_folder.name
    run_name = sample_folder.parts[-2]
    date = run_name.split("_")[0]

    html_report = read_file_as_blob(
        sample_folder / "REPORT" / f"{sample_name}.html"
    )

    fastp_json_path = sample_folder / "FASTP" / f"{sample_name}.fastp.json"
    with open(fastp_json_path) as fh:
        fastp_summary = json.load(fh).get("summary", {})
    before = fastp_summary.get("before_filtering", {})
    read_len = before.get("read1_mean_length", "")
    number_reads = before.get("total_reads", 0)

    flagstat_proc = FlagstatProcessor()
    flagstat_path = sample_folder / "logs" / "human_contamination_flagstat.txt"
    flagstat_df = flagstat_proc.process(str(flagstat_path))
    flagstat_lookup = dict(zip(flagstat_df["metric"], flagstat_df["value"]))
    percent_mapped = flagstat_lookup.get("percent_mapped", 0.0)

    # PCR duplicate stats from `samtools markdup -s` on the host BAM.
    # The file is optional (missing on legacy runs); when absent both
    # columns are left blank so a parity diff against pre-markdup runs
    # is still clean once the new columns are dropped.
    duplicate_pairs: int | str = ""
    duplicate_rate_percent: float | str = ""
    markdup_path = sample_folder / "logs" / "human_markdup_stats.txt"
    if markdup_path.exists():
        markdup_stats: dict[str, int] = {}
        with open(markdup_path) as fh:
            for line in fh:
                if ":" not in line:
                    continue
                key, _, value = line.strip().partition(":")
                try:
                    markdup_stats[key.strip()] = int(value.strip())
                except ValueError:
                    # Skip non-integer values (e.g. the COMMAND header line).
                    continue
        examined = markdup_stats.get("EXAMINED", 0)
        total_dup = markdup_stats.get("DUPLICATE TOTAL", 0)
        duplicate_pairs = markdup_stats.get("DUPLICATE PAIR", 0)
        duplicate_rate_percent = (
            100.0 * total_dup / examined if examined > 0 else 0.0
        )

    # Kraken Domain-level viral percent. The Kraken wrangled CSV
    # contains a single row with name == "Viruses"; its percent
    # already accounts for every species clade beneath it.
    # tax_lvl is "D" in the standard pluspf DB and "R1" in the small
    # viral-only DBs (e.g. k2_viral_*), so accept either.
    kraken_df = pd.read_csv(
        sample_folder / "KRAKEN" / f"{sample_name}.kraken.csv"
    )
    domain_rows = kraken_df.loc[
        (kraken_df["tax_lvl"].isin(["D", "R1"])) & (kraken_df["name"] == "Viruses"),
        "percent",
    ]
    kraken_virus_percent = float(domain_rows.iloc[0]) if not domain_rows.empty else 0.0

    # Kaiju: drop unclassified / "cannot be assigned" rows (taxon_id NA in
    # the table) before summing percents, matching the original
    # virusHanter behaviour. Keep them in the table used for the
    # top-N name list, which the original also included.
    kaiju_table_path = sample_folder / "KAIJU" / f"{sample_name}.kaiju.table.tsv"
    kaiju_report = read_file_as_blob(kaiju_table_path)
    kaiju_df = pd.read_csv(kaiju_table_path, sep="\t")
    kaiju_virus_percent = float(kaiju_df.dropna()["percent"].sum())
    top_virus_kaiju = "||".join(
        f"{row.taxon_name} ({row.reads})"
        for row in kaiju_df.head(10).itertuples()
    )

    blastn_csv = sample_folder / "BLASTN" / f"{sample_name}.contigs.blastn.csv"
    if blastn_csv.exists():
        blastn_df = pd.read_csv(blastn_csv)
        blastn_report = read_file_as_blob(blastn_csv)
    else:
        blastn_df = pd.DataFrame()
        blastn_report = ""

    number_contigs = len(blastn_df)
    if {"match_name", "read_len"}.issubset(blastn_df.columns):
        top_contigs_blastn = "||".join(
            f"{row.match_name} ({row.read_len})"
            for row in blastn_df.head(5).itertuples()
        )
    elif "match_name" in blastn_df.columns:
        top_contigs_blastn = "||".join(
            blastn_df["match_name"].head(5).astype(str).tolist()
        )
    else:
        top_contigs_blastn = ""

    return pd.DataFrame([{
        "run_name": run_name,
        "sample_name": sample_name,
        "date": date,
        "read_len": read_len,
        "number_reads": number_reads,
        "mapped_to_human_percent": percent_mapped,
        "kraken_virus_percent": kraken_virus_percent,
        "kaiju_virus_percent": kaiju_virus_percent,
        "number_of_contigs": number_contigs,
        "top_contigs_blastn": top_contigs_blastn,
        "top_virus_kaiju": top_virus_kaiju,
        "html_report": html_report,
        "kaiju_report": kaiju_report,
        "blastn_report": blastn_report,
        # Trailing columns added 2026-05-17 (Twist VRP audit). Left blank
        # on legacy runs that pre-date the markdup_human rule; a
        # column-dropped diff against an older run should still be clean.
        "duplicate_pairs": duplicate_pairs,
        "duplicate_rate_percent": duplicate_rate_percent,
    }])


def main() -> None:
    # snakemake is injected into globals by the Snakemake `script:` runner.
    sm = globals()["snakemake"]  # noqa: F821 (provided at runtime)
    results_folder = Path(sm.params.results_folder)
    samples = [Path(report).parent.parent.name for report in sm.input.reports]

    rows = [aggregate_sample_info(results_folder / sample) for sample in samples]
    run_info_df = pd.concat(rows, ignore_index=True)
    run_info_df.to_csv(sm.output.run_info_csv, index=False)


if __name__ == "__main__":
    main()
