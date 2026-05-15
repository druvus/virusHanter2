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
    sample_name = sample_folder.name

    report_html = read_file_as_blob(
        sample_folder / "REPORT" / f"{sample_name}.html"
    )

    fastp_json_path = sample_folder / "FASTP" / f"{sample_name}.fastp.json"
    with open(fastp_json_path) as fh:
        fastp_summary = json.load(fh).get("summary", {})
    before = fastp_summary.get("before_filtering", {})
    read_length = before.get("read1_mean_length", "")
    number_reads = before.get("total_reads", 0)

    flagstat_proc = FlagstatProcessor()
    flagstat_path = sample_folder / "logs" / "human_contamination_flagstat.txt"
    flagstat_df = flagstat_proc.process(str(flagstat_path))
    flagstat_lookup = dict(zip(flagstat_df["metric"], flagstat_df["value"]))
    percent_mapped = flagstat_lookup.get("percent_mapped", 0.0)

    kraken_df = pd.read_csv(
        sample_folder / "KRAKEN" / f"{sample_name}.kraken.csv"
    )
    kraken_virus_percent = kraken_df.loc[
        kraken_df["domain"] == "Viruses", "percent"
    ].sum()

    kaiju_df = pd.read_csv(
        sample_folder / "KAIJU" / f"{sample_name}.kaiju.table.tsv",
        sep="\t",
    )
    kaiju_virus_percent = float(kaiju_df["percent"].sum())
    top_virus_kaiju = "||".join(
        kaiju_df["taxon_name"].head(10).astype(str).tolist()
    )

    blastn_csv = sample_folder / "BLASTN" / f"{sample_name}.contigs.blastn.csv"
    blastn_df = pd.read_csv(blastn_csv) if blastn_csv.exists() else pd.DataFrame()
    number_contigs = len(blastn_df)
    top_contigs_blastn = (
        "||".join(blastn_df["match_name"].head(5).astype(str).tolist())
        if "match_name" in blastn_df.columns
        else ""
    )

    return pd.DataFrame([{
        "sample_name": sample_name,
        "read_length": read_length,
        "number_reads": number_reads,
        "mapped_to_human_percent": percent_mapped,
        "kraken_virus_percent": kraken_virus_percent,
        "kaiju_virus_percent": kaiju_virus_percent,
        "number_of_contigs": number_contigs,
        "top_contigs_blastn": top_contigs_blastn,
        "top_virus_kaiju": top_virus_kaiju,
        "report_html_blob": report_html,
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
