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

import reporthanter
from reporthanter import FlagstatProcessor

# Pull in the pipeline-side helper for hex-encoding the HTML blob.
from scripts.functions import read_file_as_blob

# Reference-database build identity (robust stamp over mtime) and
# short-path rendering live in scripts/provenance.py so the sidecar
# writer reuses exactly the same logic. Tool versions come from the
# resolved-version collector.
from scripts.collect_software_versions import headline_versions
from scripts.provenance import (
    databases_build_identity_string,
    databases_provenance_span_days,
    databases_provenance_string,
    databases_used_string,
    db_build_identity,
)


def _tool_versions_string(software_versions_tsv: Path | None) -> str:
    """Render the resolved headline tool versions as a compact
    ``fastp=0.24.0;bwa=0.7.18;...`` cell. Empty when the collector
    output is missing (e.g. a legacy run) so a column-dropped parity
    diff stays clean.
    """
    if software_versions_tsv is None or not Path(software_versions_tsv).is_file():
        return ""
    rows: list[dict[str, str]] = []
    with Path(software_versions_tsv).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            if len(values) == len(header):
                rows.append(dict(zip(header, values)))
    headline = headline_versions(rows)
    return ";".join(f"{tool}={ver}" for tool, ver in sorted(headline.items()))


def _parse_quast_report(report_tsv: Path) -> dict[str, float | int]:
    """Pull n_contigs and N50 out of a QUAST ``report.tsv``.

    QUAST writes a two-column TSV: ``Assembly`` then a single sample
    column. Returns an empty dict on a missing or malformed file so
    the caller can degrade gracefully.
    """
    out: dict[str, float | int] = {}
    if not report_tsv.exists() or report_tsv.stat().st_size == 0:
        return out
    try:
        with open(report_tsv) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                key, value = parts[0].strip(), parts[1].strip()
                if key == "# contigs":
                    try:
                        out["n_contigs"] = int(value)
                    except ValueError:
                        continue
                elif key == "N50":
                    try:
                        out["n50"] = int(value)
                    except ValueError:
                        continue
    except OSError:
        return {}
    return out


def aggregate_sample_info(
    sample_folder: Path,
    *,
    databases_used: str = "",
    databases_provenance: str = "",
    databases_span_days: int = 0,
    databases_build_identity: str = "",
    tool_versions: str = "",
    reporthanter_version: str = "",
    assemblers: list[str] | None = None,
) -> pd.DataFrame:
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
    # The bwa backend writes ``human_contamination_flagstat.txt``;
    # the hostile backend writes ``hostile_contamination_flagstat.txt``
    # in samtools-flagstat shape. Pick whichever exists so the
    # column values stay identical regardless of host-removal tool.
    bwa_flagstat = sample_folder / "logs" / "human_contamination_flagstat.txt"
    hostile_flagstat = sample_folder / "logs" / "hostile_contamination_flagstat.txt"
    if hostile_flagstat.exists():
        flagstat_path = hostile_flagstat
        host_removal_tool = "hostile"
    else:
        flagstat_path = bwa_flagstat
        host_removal_tool = "bwa"
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

    # Read every per-assembler BLASTN CSV under
    # `{sample_folder}/{assembler}/BLASTN/...` and concatenate. The
    # `assembler` column comes through from `wrangle_pilon`. A blob of
    # the concatenated CSV is stored under `blastn_report` so the
    # parity-locked column still carries the full BLAST table.
    asm_list = list(assemblers) if assemblers else ["MEGAHIT"]
    blastn_frames: list[pd.DataFrame] = []
    for asm in asm_list:
        csv = sample_folder / asm / "BLASTN" / f"{sample_name}.contigs.blastn.csv"
        if csv.exists() and csv.stat().st_size > 0:
            frame = pd.read_csv(csv)
            if "assembler" not in frame.columns:
                frame = frame.assign(assembler=asm)
            blastn_frames.append(frame)
    blastn_df = (
        pd.concat(blastn_frames, ignore_index=True)
        if blastn_frames
        else pd.DataFrame()
    )

    # `blastn_report` keeps the parity-locked hex-encoded blob shape.
    # When multiple assemblers contribute, persist the concatenated
    # CSV to a temp path long enough to read_file_as_blob it.
    if not blastn_df.empty:
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as tmp:
            blastn_df.to_csv(tmp.name, index=False)
            blastn_report = read_file_as_blob(Path(tmp.name))
            Path(tmp.name).unlink(missing_ok=True)
    else:
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

    # Per-assembler trailing diagnostics: contig count, QUAST n_contigs
    # / N50 (when QUAST is enabled). Always emitted so the column set
    # stays stable across samples within a batch, even when a sample
    # produced zero contigs from one assembler.
    per_asm_columns: dict[str, int | float | str] = {}
    for asm in asm_list:
        if not blastn_df.empty and "assembler" in blastn_df.columns:
            per_asm_columns[f"{asm.lower()}_n_contigs"] = int(
                (blastn_df["assembler"] == asm).sum()
            )
        else:
            per_asm_columns[f"{asm.lower()}_n_contigs"] = 0
        quast_tsv = sample_folder / asm / "QUAST" / "report.tsv"
        quast_stats = _parse_quast_report(quast_tsv)
        per_asm_columns[f"{asm.lower()}_n50"] = quast_stats.get("n50", "")
    assemblers_used = ";".join(asm_list)

    # Optional geNomad summary. Only present when the workflow ran
    # with GENOMAD: "TRUE". Sums across per-assembler summaries when
    # multi-assembler mode is on; the column shape is unchanged so a
    # column-dropped diff against pre-multi-assembler runs is clean.
    genomad_viral_contigs: int | str = ""
    genomad_max_virus_score: float | str = ""
    genomad_totals: list[int] = []
    genomad_maxes: list[float] = []
    for asm in asm_list:
        # geNomad names its outputs after the input FASTA stem
        # (<sample>_improved_contigs.fasta), so the rule writes
        # <asm>/GENOMAD/<sample>_improved_contigs_summary/<sample>_improved_contigs_virus_summary.tsv
        # (matching rules/assembly.smk:480 and the per_virus_metrics /
        # generate_report consumers in post_processing.smk). The earlier
        # <sample>_summary/<sample>_virus_summary.tsv form never resolved,
        # leaving these columns silently blank on GENOMAD: "TRUE" runs.
        genomad_path = (
            sample_folder
            / asm
            / "GENOMAD"
            / f"{sample_name}_improved_contigs_summary"
            / f"{sample_name}_improved_contigs_virus_summary.tsv"
        )
        if genomad_path.exists() and genomad_path.stat().st_size > 0:
            try:
                gdf = pd.read_csv(genomad_path, sep="\t")
                if "virus_score" in gdf.columns:
                    genomad_totals.append(int(len(gdf)))
                    if len(gdf):
                        genomad_maxes.append(float(gdf["virus_score"].max()))
            except Exception:  # noqa: BLE001
                continue
    if genomad_totals:
        genomad_viral_contigs = int(sum(genomad_totals))
    if genomad_maxes:
        genomad_max_virus_score = float(max(genomad_maxes))

    row = {
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
        # Trailing columns added 2026-05-21 (geNomad surface). Only
        # populated when GENOMAD: "TRUE"; otherwise empty so a
        # column-dropped diff against pre-geNomad runs is clean.
        "genomad_viral_contigs": genomad_viral_contigs,
        "genomad_max_virus_score": genomad_max_virus_score,
        # Trailing diagnostic-provenance columns. `databases_used`
        # records the reference DB paths the workflow consumed (one
        # KEY=VALUE per DB, semicolon-separated). `reporthanter_version`
        # captures the version that rendered the per-sample HTML
        # blobs. Both columns are constant within a batch but vary
        # across batches and across DB / reporthanter refreshes.
        "databases_used": databases_used,
        # `databases_provenance` carries one `KEY=YYYY-MM-DD` per
        # configured DB (from each DB's representative-file mtime
        # or — for VIRUS_PARQUET — the explicit `build_date_utc` in
        # its build_stats.json sidecar). `databases_span_days` is
        # the span between the oldest and newest of those dates;
        # >180 means classifier DBs likely came from divergent
        # snapshots and the cross-DB consistency the refresh
        # workflow assumes is no longer guaranteed.
        "databases_provenance": databases_provenance,
        "databases_span_days": databases_span_days,
        # `databases_build_identity` carries one `KEY=<identity>` per DB,
        # preferring a robust build stamp (the refresh workflow's
        # build_stats.json source+date, or the version-bearing directory
        # name such as `checkv-db-v1.5` / `k2_pluspf_20240112`) over a
        # bare mtime. `tool_versions` records the conda-resolved version
        # of each headline tool that actually ran
        # (`fastp=0.24.0;bwa=0.7.18;...`). Both are blank on legacy runs
        # so a column-dropped parity diff stays clean.
        "databases_build_identity": databases_build_identity,
        "tool_versions": tool_versions,
        "reporthanter_version": reporthanter_version,
        # Trailing multi-assembler diagnostics. `assemblers_used` is a
        # semicolon-delimited list of assemblers run for this sample;
        # per-assembler trailing columns (`{assembler}_n_contigs`,
        # `{assembler}_n50`) carry the equivalent of QUAST's headline
        # numbers when QUAST is enabled. Empty / zero on runs that did
        # not enable QUAST for that assembler.
        "assemblers_used": assemblers_used,
        # Trailing column recording which host-removal tool produced
        # the per-sample read set ("bwa" or "hostile"). Allows
        # downstream comparisons of duplicate / host-carryover rates
        # across runs that switched backends.
        "host_removal_tool": host_removal_tool,
    }
    row.update(per_asm_columns)
    return pd.DataFrame([row])


def main() -> None:
    # snakemake is injected into globals by the Snakemake `script:` runner.
    sm = globals()["snakemake"]  # noqa: F821 (provided at runtime)
    results_folder = Path(sm.params.results_folder)
    samples = [Path(report).parent.parent.name for report in sm.input.reports]

    cfg = dict(sm.config)
    identity = db_build_identity(cfg)
    databases_used = databases_used_string(identity)
    databases_provenance = databases_provenance_string(identity)
    databases_span_days = databases_provenance_span_days(identity)
    databases_build_identity = databases_build_identity_string(identity)
    if databases_span_days > 180:
        print(
            "[aggregate_run_information] WARNING: reference DB build "
            f"dates span {databases_span_days} days "
            f"({databases_provenance}); cross-DB taxonomy may be out "
            "of sync. Re-run refresh/refresh_virus_parquet.smk to "
            "rebuild VIRUS_PARQUET and KAIJU_DB from a current "
            "snapshot."
        )
    software_versions_tsv = getattr(sm.input, "software_versions", None)
    tool_versions = _tool_versions_string(
        Path(software_versions_tsv) if software_versions_tsv else None
    )
    reporthanter_version = getattr(reporthanter, "__version__", "")
    assemblers = list(sm.params.get("assemblers", ["MEGAHIT"]))

    rows = [
        aggregate_sample_info(
            results_folder / sample,
            databases_used=databases_used,
            databases_provenance=databases_provenance,
            databases_span_days=databases_span_days,
            databases_build_identity=databases_build_identity,
            tool_versions=tool_versions,
            reporthanter_version=reporthanter_version,
            assemblers=assemblers,
        )
        for sample in samples
    ]
    run_info_df = pd.concat(rows, ignore_index=True)
    run_info_df.to_csv(sm.output.run_info_csv, index=False)


if __name__ == "__main__":
    main()
