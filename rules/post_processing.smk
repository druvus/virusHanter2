# post_processing.smk
#
# Coverage alignment, report generation (via the reporthanter CLI), run
# aggregation, and optional cleanup.

import re

from scripts.functions import read_file_as_blob

# Set variables
THREADS = config["THREADS"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(config["SAMPLES"]).name)
VIRUS_PARQUET = config["VIRUS_PARQUET"]
NUMBER_OF_PLOTS = config["NUMBER_OF_PLOTS"]
# Window size (bp) for mosdepth's --by flag. Smaller values give a
# higher-resolution coverage trace in the report at the cost of a
# larger regions.bed.gz; 100 bp is a sensible default for the viral
# references this pipeline targets.
COVERAGE_WINDOW = int(config.get("COVERAGE_WINDOW", 100))
SECONDARY_HOST_NAME = config.get("SECONDARY_HOST_NAME", "")

# Rule: Build a per-sample reference set from the union of every
# enabled classifier's viral hits and align reads to it.
#
# Sources (each can be enabled or disabled via config[COVERAGE_SOURCES]):
#   - KRAKEN: top-N viral taxa from the wrangled Kraken2 CSV.
#   - KAIJU: top-N viral taxa from the kaiju2table output, filtered
#     against VIRUS_PARQUET so non-viral RefSeq hits (Kaiju's default
#     refseq DB is broader than viral) do not enter.
#   - BLAST: every per-assembler merged BLAST CSV; taxids are
#     resolved through VIRUS_PARQUET's accession lookup.
#
# Outputs include a `unmapped_taxids.tsv` sidecar listing classified
# taxids that VIRUS_PARQUET has no reference for, so the reviewer can
# see why a hit did not produce a coverage trace. The `virus_names`
# sidecar gains a `sources` column tagging each chrom with the
# classifier(s) that contributed it.
rule bwa_align_to_kraken_hits:
    input:
        kraken_csv=rules.wrangle_kraken.output.kraken_csv,
        kaiju_table=rules.kaiju_to_table.output.kaiju_table,
        blastn_csvs=expand(
            f"{RESULT_FOLDER}/{{{{sample}}}}/{{assembler}}/CHECKV/{{{{sample}}}}.merged.csv",
            assembler=ASSEMBLERS,
        ),
        r1=lambda wildcards: rules.bam_to_fastq_human.output.r1 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r1,
        r2=lambda wildcards: rules.bam_to_fastq_human.output.r2 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r2,
    output:
        virus_fasta=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/kraken_top_viruses.fasta",
        virus_names=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/kraken_top_virus_names.tsv",
        unmapped_taxids=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/unmapped_taxids.tsv",
        bam=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/{{sample}}_kraken.bam",
    params:
        virus_db=VIRUS_PARQUET,
        index_prefix=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/bwa/{{sample}}",
        coverage_sources=COVERAGE_SOURCES,
        coverage_top_n=COVERAGE_TOP_N,
        taxdump_nodes=TAXDUMP_NODES,
        coverage_rank_filter=COVERAGE_RANK_FILTER,
        coverage_genus_walkup=COVERAGE_GENUS_WALKUP,
    threads: THREADS
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_kraken.log"
    conda:
        "../envs/bwa.yaml"
    run:
        import pandas as pd
        from pathlib import Path

        from scripts.functions import parquet_accession_to_taxid
        from scripts.build_virus_parquet import (
            find_genus_taxid,
            parse_nodes_dmp,
        )

        # Load the reference parquet once. Subsequent lookups are
        # by tax_id and by accession.
        virus_db_df = pd.read_parquet(params.virus_db)
        parquet_tax_ids = set(
            virus_db_df["tax_id"].dropna().astype(int).tolist()
        )
        acc_to_tax = parquet_accession_to_taxid(virus_db_df)

        # Inverse index from genus_taxid to a representative parquet
        # tax_id under that genus. Used by the walk-up fallback when
        # a classifier hit is absent from the parquet AND the genus
        # itself is also absent (typical for viral RefSeq, which
        # has no genus-level references — only species-level rows
        # that carry the genus as a column). We pick the longest
        # sequence under each genus as the representative.
        genus_to_rep_taxid: dict[int, int] = {}
        if "genus_taxid" in virus_db_df.columns:
            sub = virus_db_df.loc[virus_db_df["genus_taxid"] > 0].copy()
            if not sub.empty:
                sub = sub.assign(_seqlen=sub["sequence"].str.len())
                sub = sub.sort_values(
                    ["genus_taxid", "_seqlen"],
                    ascending=[True, False],
                    kind="mergesort",
                )
                for r in sub.drop_duplicates(subset=["genus_taxid"], keep="first").itertuples():
                    genus_to_rep_taxid[int(r.genus_taxid)] = int(r.tax_id)

        # Optional taxdump-driven rank lookup. When TAXDUMP_NODES is
        # empty or the file is missing, the rank filter and walk-up
        # both degrade to no-ops; the rule still produces the bare
        # multi-source union as before.
        rank_filter: set[str] = set(params.coverage_rank_filter or [])
        taxdump_path = params.taxdump_nodes
        nodes: dict[int, tuple[int, str]] = {}
        if taxdump_path and Path(taxdump_path).is_file():
            nodes = parse_nodes_dmp(Path(taxdump_path))
        elif rank_filter or params.coverage_genus_walkup:
            print(
                "[bwa_align_to_kraken_hits] TAXDUMP_NODES not set or missing; "
                "rank filter and genus walk-up are disabled."
            )
            rank_filter = set()

        # Per-source taxid sets, plus a parallel name map keyed by
        # tax_id so the eventual sidecar can carry the species name
        # the classifier originally reported.
        sources_for_tid: dict[int, set[str]] = {}
        names_for_tid: dict[int, str] = {}
        unmapped_rows: list[tuple[int, str, str, str]] = []  # (tid, name, source, reason)

        def _record(tid: int, name: str, source: str) -> None:
            rank = nodes.get(tid, (0, "unknown"))[1] if nodes else "unknown"
            # Higher-rank propagation rows (kingdom, phylum, class,
            # order, family) are dropped silently — they have no
            # per-taxid reference and would otherwise flood the
            # unmapped sidecar.
            if rank in rank_filter:
                return
            if tid in parquet_tax_ids:
                sources_for_tid.setdefault(tid, set()).add(source)
                # Prefer the first non-empty name we see; classifiers
                # report the same species under cosmetically different
                # strings. KRAKEN is processed first so its naming wins.
                if name and tid not in names_for_tid:
                    names_for_tid[tid] = name
                return
            # Walk up to the genus and substitute a representative
            # parquet reference for that genus. Two-step:
            #   1. Resolve the genus_tid via the taxdump.
            #   2. Look the genus_tid up in the parquet by tax_id
            #      first, then by ``genus_taxid`` column (the typical
            #      case for viral RefSeq, which has no genus-level
            #      references — only species rows tagged with their
            #      genus).
            # The substitution tags the source with ``->genus`` so
            # the reviewer sees the fallback in the coverage tab
            # label.
            if params.coverage_genus_walkup and nodes:
                genus_tid = find_genus_taxid(tid, nodes)
                if genus_tid:
                    if genus_tid in parquet_tax_ids:
                        rep = genus_tid
                    else:
                        rep = genus_to_rep_taxid.get(genus_tid, 0)
                    if rep:
                        sources_for_tid.setdefault(rep, set()).add(
                            f"{source}->genus"
                        )
                        if name and rep not in names_for_tid:
                            names_for_tid[rep] = name
                        return
            unmapped_rows.append((tid, name, source, "absent_from_parquet"))

        if "KRAKEN" in params.coverage_sources:
            kraken_df = pd.read_csv(input.kraken_csv)
            top_kraken = (
                kraken_df.loc[kraken_df.domain == "Viruses"]
                .sort_values("percent", ascending=False)
                .head(int(params.coverage_top_n))
            )
            for r in top_kraken.itertuples():
                try:
                    tid = int(r.taxonomy_id)
                except (ValueError, TypeError):
                    continue
                _record(tid, str(getattr(r, "name", "")), "kraken")

        if "KAIJU" in params.coverage_sources:
            try:
                kaiju_df = pd.read_csv(input.kaiju_table, sep="\t")
            except Exception:  # noqa: BLE001
                kaiju_df = pd.DataFrame()
            if not kaiju_df.empty and "taxon_id" in kaiju_df.columns:
                kaiju_df = kaiju_df.dropna(subset=["taxon_id"])
                kaiju_df = kaiju_df.loc[
                    kaiju_df["taxon_name"].fillna("") != "unclassified"
                ]
                if "percent" in kaiju_df.columns:
                    kaiju_df = kaiju_df.sort_values("percent", ascending=False)
                for r in kaiju_df.head(int(params.coverage_top_n)).itertuples():
                    try:
                        tid = int(r.taxon_id)
                    except (ValueError, TypeError):
                        continue
                    _record(tid, str(getattr(r, "taxon_name", "")), "kaiju")

        if "BLAST" in params.coverage_sources:
            for csv in input.blastn_csvs:
                if not Path(csv).exists() or Path(csv).stat().st_size == 0:
                    continue
                try:
                    blast_df = pd.read_csv(csv)
                except Exception:  # noqa: BLE001
                    continue
                if "accession" not in blast_df.columns:
                    continue
                seen_per_csv: set[int] = set()
                for r in blast_df.itertuples():
                    accession = getattr(r, "accession", None)
                    if accession is None or pd.isna(accession):
                        continue
                    acc_str = str(accession).strip()
                    tid = acc_to_tax.get(acc_str) or acc_to_tax.get(
                        acc_str.split(".")[0]
                    )
                    if tid is None:
                        continue
                    if tid in seen_per_csv:
                        continue
                    seen_per_csv.add(tid)
                    _record(tid, str(getattr(r, "match_name", "")), "blast")

        # Resolve the parquet rows for every retained tax_id.
        selected_viruses = virus_db_df[
            virus_db_df["tax_id"].astype(int).isin(sources_for_tid.keys())
        ]

        Path(output.virus_fasta).parent.mkdir(parents=True, exist_ok=True)
        with open(output.virus_fasta, "w") as f, open(output.virus_names, "w") as nf:
            nf.write("chrom\ttax_id\tname\tsources\n")
            for row in selected_viruses.itertuples():
                accession = row.name.strip().split()[0]
                tid = int(row.tax_id)
                species = names_for_tid.get(tid, "")
                source_tag = ";".join(sorted(sources_for_tid.get(tid, set())))
                f.write(f">{row.name.strip()}\n{row.sequence}\n")
                nf.write(f"{accession}\t{tid}\t{species}\t{source_tag}\n")

        with open(output.unmapped_taxids, "w") as uf:
            uf.write("tax_id\tname\tsource\treason\n")
            for tid, name, source, reason in unmapped_rows:
                # Strip tab/newline from name to keep the TSV well-formed.
                clean_name = name.replace("\t", " ").replace("\n", " ")
                uf.write(f"{tid}\t{clean_name}\t{source}\t{reason}\n")

        # If no taxid landed (a smoke-test scenario with empty
        # classifiers), emit a single dummy reference so BWA still
        # produces a valid index + BAM and downstream rules do not
        # crash on an empty FASTA.
        if Path(output.virus_fasta).stat().st_size == 0:
            with open(output.virus_fasta, "w") as f:
                f.write(">DUMMY_REF\n")
                f.write("N" * 100 + "\n")

        # Create BWA index
        index_prefix = params.index_prefix
        Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)
        shell("bwa index -p {index_prefix} {output.virus_fasta} > {log} 2>&1")

        # Align reads to viral sequences
        shell("bwa mem -t {threads} {index_prefix} {input.r1} {input.r2} | samtools sort -o {output.bam} - >> {log} 2>&1")
        shell("samtools index {output.bam}")

        # Clean up index files
        shell("rm -rf {index_prefix}*")

# Rule: Per-reference coverage statistics from the kraken-top viral BAM.
# Output is a mosdepth regions BED keyed by reference, consumed by
# reporthanter to render interactive Altair coverage traces in the
# per-sample HTML report.
rule mosdepth_kraken_hits:
    input:
        bam=rules.bwa_align_to_kraken_hits.output.bam,
    output:
        summary=f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}.mosdepth.summary.txt",
        regions=f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}.regions.bed.gz",
        # `--thresholds 1,5,10` emits a sibling `thresholds.bed.gz`
        # with per-region counts of bases at each coverage threshold.
        # `per_virus_metrics` sums the 5x column per chrom to get
        # `bases_above_5x`.
        thresholds=f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}.thresholds.bed.gz",
    params:
        prefix=f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}",
        window=COVERAGE_WINDOW,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/mosdepth.log"
    threads: 4
    conda:
        "../envs/mosdepth.yaml"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        mosdepth -t {threads} --no-per-base --by {params.window} \
            --thresholds 1,5,10 \
            {params.prefix} {input.bam} > {log} 2>&1
        """


# Rule: Generate interactive report via the reporthanter CLI
#
# Consumes the per-(sample, assembler) merged CSVs and any optional
# per-assembler QUAST / geNomad outputs and hands them to the
# reporthanter CLI. The CLI accepts repeated --blastn_file,
# --quast_report and --genomad_summary flags so the report carries an
# `assembler` column on the contig table.
rule generate_report:
    input:
        flagstat=rules.remove_host.output.flagstat,
        secondary_flagstat=rules.remove_secondary_host.output.flagstat,
        fastp_json=rules.fastp.output.json_report,
        blastn_csvs=expand(
            f"{RESULT_FOLDER}/{{{{sample}}}}/{{assembler}}/CHECKV/{{{{sample}}}}.merged.csv",
            assembler=ASSEMBLERS,
        ),
        # reporthanter's KrakenProcessor reads the raw Kraken2 report
        # (TSV, 6 columns, no header) and wrangles internally. The
        # pipeline-side wrangle_kraken CSV is used by
        # aggregate_run_information instead.
        kraken_report=rules.kraken.output.kraken_report,
        kaiju_table=rules.kaiju_to_table.output.kaiju_table,
        mosdepth_regions=rules.mosdepth_kraken_hits.output.regions,
        virus_names=rules.bwa_align_to_kraken_hits.output.virus_names,
        # When QUAST is enabled, surface every per-assembler report
        # inside the HTML as an Alignment Stats sub-tab.
        **(
            {
                "quast_reports": expand(
                    f"{RESULT_FOLDER}/{{{{sample}}}}/{{assembler}}/QUAST/report.tsv",
                    assembler=ASSEMBLERS,
                )
            }
            if config.get("QUAST", "FALSE") == "TRUE"
            else {}
        ),
        # When geNomad is enabled, surface every per-assembler virus
        # summary as a Classification of Contigs sub-tab.
        **(
            {
                "genomad_summaries": expand(
                    f"{RESULT_FOLDER}/{{{{sample}}}}/{{assembler}}/GENOMAD/{{{{sample}}}}_summary/{{{{sample}}}}_virus_summary.tsv",
                    assembler=ASSEMBLERS,
                )
            }
            if config.get("GENOMAD", "FALSE") == "TRUE"
            else {}
        ),
    output:
        report_html=f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html",
    conda:
        "../envs/reporthanter.yaml"
    params:
        secondary_args=(
            lambda wildcards: (
                f"--secondary_flagstat_file {RESULT_FOLDER}/{wildcards.sample}/logs/secondary_contamination_flagstat.txt "
                f"--secondary_host {SECONDARY_HOST_NAME}"
            )
            if SECONDARY_HOST_OR_NOT
            else ""
        ),
        # Strip the trailing "_R" that the paired-read wildcard scheme
        # leaves behind (Illumina-style R1/R2 file pairs split into a
        # sample name ending in "_R"). The wildcard itself must stay as
        # is for filesystem paths; this affects only the display name
        # shown in the report header.
        display_name=lambda wildcards: re.sub(r"_R$", "", wildcards.sample),
        blastn_args=lambda wildcards, input: " ".join(
            f"--blastn_file {p}" for p in input.blastn_csvs
        ),
        quast_args=lambda wildcards, input: (
            " ".join(f"--quast_report {p}" for p in input.quast_reports)
            if hasattr(input, "quast_reports")
            else ""
        ),
        genomad_args=lambda wildcards, input: (
            " ".join(
                f"--genomad_summary {p}" for p in input.genomad_summaries
            )
            if hasattr(input, "genomad_summaries")
            else ""
        ),
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/reporthanter.log",
    shell:
        """
        reporthanter \
            {params.blastn_args} \
            --kraken_file {input.kraken_report} \
            --kaiju_table {input.kaiju_table} \
            --fastp_json {input.fastp_json} \
            --flagstat_file {input.flagstat} \
            --mosdepth_regions {input.mosdepth_regions} \
            --virus_names {input.virus_names} \
            {params.quast_args} \
            {params.genomad_args} \
            --output {output.report_html} \
            --sample_name {params.display_name} \
            {params.secondary_args} \
            > {log} 2>&1
        """

# Rule: Per-sample per-virus metrics.
#
# Joins the existing pipeline outputs (Kraken, Kaiju, BLASTN merged
# CSV, mosdepth summary + thresholds, fastp JSON, host flagstat) and
# the workflow-level viral parquet into a flat CSV with one row per
# detected Kraken viral taxid for this sample. Schema: see
# `docs/PER_VIRUS_OUTPUT.md`.
rule per_virus_metrics:
    input:
        kraken_csv=rules.wrangle_kraken.output.kraken_csv,
        kaiju_tsv=rules.kaiju_to_table.output.kaiju_table,
        blastn_csvs=expand(
            f"{RESULT_FOLDER}/{{{{sample}}}}/{{assembler}}/CHECKV/{{{{sample}}}}.merged.csv",
            assembler=ASSEMBLERS,
        ),
        mosdepth_summary=rules.mosdepth_kraken_hits.output.summary,
        mosdepth_thresholds=rules.mosdepth_kraken_hits.output.thresholds,
        fastp_json=rules.fastp.output.json_report,
        flagstat=rules.remove_host.output.flagstat,
        virus_parquet=VIRUS_PARQUET,
        # When geNomad is enabled, take every per-assembler summary as
        # an extra input. Evaluated at rule-build time.
        **(
            {
                "genomad_summaries": expand(
                    f"{RESULT_FOLDER}/{{{{sample}}}}/{{assembler}}/GENOMAD/{{{{sample}}}}_summary/{{{{sample}}}}_virus_summary.tsv",
                    assembler=ASSEMBLERS,
                )
            }
            if config.get("GENOMAD", "FALSE") == "TRUE"
            else {}
        ),
    output:
        per_virus_csv=f"{RESULT_FOLDER}/{{sample}}/{{sample}}.per_virus.csv",
    params:
        run_name=Path(config["SAMPLES"]).name,
        top_n=NUMBER_OF_PLOTS,
        blastn_args=lambda wildcards, input: " ".join(
            f"--blastn-csv {p}" for p in input.blastn_csvs
        ),
        genomad_args=lambda wildcards, input: (
            " ".join(
                f"--genomad-summary {p}" for p in input.genomad_summaries
            )
            if hasattr(input, "genomad_summaries")
            else ""
        ),
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/per_virus_metrics.log"
    conda:
        "../envs/panel.yaml"
    shell:
        """
        python scripts/per_virus_metrics.py \
            --sample-name {wildcards.sample} \
            --run-name {params.run_name} \
            --kraken-csv {input.kraken_csv} \
            --kaiju-tsv {input.kaiju_tsv} \
            {params.blastn_args} \
            --mosdepth-summary {input.mosdepth_summary} \
            --mosdepth-thresholds {input.mosdepth_thresholds} \
            --fastp-json {input.fastp_json} \
            --flagstat {input.flagstat} \
            --virus-parquet {input.virus_parquet} \
            --top-n {params.top_n} \
            {params.genomad_args} \
            --out {output.per_virus_csv} \
            > {log} 2>&1
        """


# Rule: Concatenate per-sample per_virus CSVs into a single batch file.
rule aggregate_per_virus:
    input:
        per_sample=expand(
            f"{RESULT_FOLDER}/{{sample}}/{{sample}}.per_virus.csv",
            sample=SAMPLES,
        ),
    output:
        per_virus_csv=f"{RESULT_FOLDER}/per_virus_{Path(config['SAMPLES']).name}.csv",
    log:
        f"{RESULT_FOLDER}/logs/aggregate_per_virus.log"
    conda:
        "../envs/panel.yaml"
    shell:
        """
        python scripts/aggregate_per_virus.py \
            --in {input.per_sample} \
            --out {output.per_virus_csv} \
            > {log} 2>&1
        """


# Rule: Workflow-level MultiQC aggregation.
# Runs after every per-sample report is finalised and after the
# aggregate CSV is written, scans RESULT_FOLDER for fastp/samtools/
# kraken/mosdepth/markdup outputs, and emits a single HTML for the
# whole batch. Gated by the MULTIQC config flag in Snakefile rule all.
rule multiqc:
    input:
        # The aggregated CSV is the latest "everything done" sentinel
        # in the workflow, so depending on it pulls in all per-sample
        # reports and stats files transitively.
        run_info_csv=f"{RESULT_FOLDER}/run_information_{Path(config['SAMPLES']).name}.csv",
        markdup=expand(f"{RESULT_FOLDER}/{{sample}}/logs/human_markdup_stats.txt", sample=SAMPLES),
        mosdepth=expand(f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}.mosdepth.summary.txt", sample=SAMPLES),
        # Optional: QUAST reports when the assembler-QC step is enabled.
        # MultiQC scans the results folder regardless, but depending on
        # the report.tsv keeps the dependency graph honest so MultiQC
        # waits until QUAST has finished before scanning. One report
        # per (sample, assembler).
        quast=(
            expand(
                f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/QUAST/report.tsv",
                sample=SAMPLES,
                assembler=ASSEMBLERS,
            )
            if config.get("QUAST", "FALSE") == "TRUE"
            else []
        ),
    output:
        html=f"{RESULT_FOLDER}/multiqc_report.html",
        data=directory(f"{RESULT_FOLDER}/multiqc_data"),
    params:
        results_folder=RESULT_FOLDER,
    log:
        f"{RESULT_FOLDER}/logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc \
            --force \
            --outdir {params.results_folder} \
            --filename multiqc_report.html \
            {params.results_folder} \
            > {log} 2>&1
        """


# Rule: Aggregate run information across samples
#
# Uses a `script:` directive (not `run:`) so the body actually executes in
# the reporthanter conda env. Snakemake's `run:` blocks always run in the
# driver Python even when a `conda:` directive is set, which would break
# the `from reporthanter import FlagstatProcessor` import.
rule aggregate_run_information:
    input:
        reports=expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
    output:
        run_info_csv=f"{RESULT_FOLDER}/run_information_{Path(config['SAMPLES']).name}.csv",
    params:
        results_folder=RESULT_FOLDER,
        assemblers=ASSEMBLERS,
    log:
        f"{RESULT_FOLDER}/logs/aggregate_run_information.log",
    conda:
        "../envs/reporthanter.yaml"
    script:
        "../scripts/aggregate_run_information.py"

# Rule: Clean up intermediate files (optional)
rule clean_everything:
    input:
        run_info_csv=rules.aggregate_run_information.output.run_info_csv,
    output:
        cleanup_done=f"{RESULT_FOLDER}/analysis_done.txt",
    params:
        results_folder=RESULT_FOLDER,
    shell:
        """
        # Remove intermediate files but keep logs, reports, CSVs, and flagstat files
        find {params.results_folder} -type f ! -name '*.html' ! -name '*.csv' ! -name '*flagstat.txt' ! -name '*.tsv' ! -name '*.log' -delete
        echo "Analysis completed on $(date)" > {output.cleanup_done}
        """