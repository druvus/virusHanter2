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

# Rule: Align reads to top Kraken2 viral hits
rule bwa_align_to_kraken_hits:
    input:
        kraken_csv=rules.wrangle_kraken.output.kraken_csv,
        r1=lambda wildcards: rules.bam_to_fastq_human.output.r1 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r1,
        r2=lambda wildcards: rules.bam_to_fastq_human.output.r2 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r2,
    output:
        virus_fasta=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/kraken_top_viruses.fasta",
        bam=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/{{sample}}_kraken.bam",
    params:
        virus_db=VIRUS_PARQUET,
        index_prefix=f"{RESULT_FOLDER}/{{sample}}/BWA_KRAKEN/bwa/{{sample}}",
    threads: THREADS
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_kraken.log"
    conda:
        "../envs/bwa.yaml"
    run:
        import pandas as pd
        from pathlib import Path

        # Read top Kraken2 viral taxa: select the 20 with the highest
        # percent classified, matching the original virusHanter ordering.
        kraken_df = pd.read_csv(input.kraken_csv)
        top_tax_ids = (
            kraken_df.loc[kraken_df.domain == "Viruses"]
            .sort_values("percent", ascending=False)
            .head(20)["taxonomy_id"]
            .tolist()
        )

        # Load viral sequences from the Parquet database
        virus_db_df = pd.read_parquet(params.virus_db)
        selected_viruses = virus_db_df[virus_db_df["tax_id"].isin(top_tax_ids)]

        # Write selected viral sequences to FASTA file
        with open(output.virus_fasta, "w") as f:
            for row in selected_viruses.itertuples():
                f.write(f">{row.name.strip()}\n{row.sequence}\n")

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
rule generate_report:
    input:
        flagstat=rules.remove_host.output.flagstat,
        secondary_flagstat=rules.remove_secondary_host.output.flagstat,
        fastp_json=rules.fastp.output.json_report,
        blastn_csv=rules.merge_checkv_blastn.output.merged_csv,
        # reporthanter's KrakenProcessor reads the raw Kraken2 report
        # (TSV, 6 columns, no header) and wrangles internally. The
        # pipeline-side wrangle_kraken CSV is used by
        # aggregate_run_information instead.
        kraken_report=rules.kraken.output.kraken_report,
        kaiju_table=rules.kaiju_to_table.output.kaiju_table,
        mosdepth_regions=rules.mosdepth_kraken_hits.output.regions,
        # When QUAST is enabled, surface its report inside the HTML
        # report as an Alignment Stats sub-tab.
        **(
            {"quast_report": rules.quast_megahit.output.report_tsv}
            if config.get("QUAST", "FALSE") == "TRUE"
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
        quast_args=(
            lambda wildcards, input: (
                f"--quast_report {input.quast_report}"
                if hasattr(input, "quast_report")
                else ""
            )
        ),
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/reporthanter.log",
    shell:
        """
        reporthanter \
            --blastn_file {input.blastn_csv} \
            --kraken_file {input.kraken_report} \
            --kaiju_table {input.kaiju_table} \
            --fastp_json {input.fastp_json} \
            --flagstat_file {input.flagstat} \
            --mosdepth_regions {input.mosdepth_regions} \
            {params.quast_args} \
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
        blastn_csv=rules.merge_checkv_blastn.output.merged_csv,
        mosdepth_summary=rules.mosdepth_kraken_hits.output.summary,
        mosdepth_thresholds=rules.mosdepth_kraken_hits.output.thresholds,
        fastp_json=rules.fastp.output.json_report,
        flagstat=rules.remove_host.output.flagstat,
        virus_parquet=VIRUS_PARQUET,
        # When geNomad is enabled, take its per-sample summary as an
        # extra input. The flag is read from `config` so the input list
        # is evaluated at rule-build time, not via a `lambda`.
        **(
            {"genomad_summary": rules.genomad.output.summary}
            if config.get("GENOMAD", "FALSE") == "TRUE"
            else {}
        ),
    output:
        per_virus_csv=f"{RESULT_FOLDER}/{{sample}}/{{sample}}.per_virus.csv",
    params:
        run_name=Path(config["SAMPLES"]).name,
        top_n=NUMBER_OF_PLOTS,
        genomad_args=(
            lambda wildcards, input: (
                f"--genomad-summary {input.genomad_summary}"
                if hasattr(input, "genomad_summary")
                else ""
            )
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
            --blastn-csv {input.blastn_csv} \
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
        # waits until QUAST has finished before scanning.
        quast=(
            expand(f"{RESULT_FOLDER}/{{sample}}/QUAST/report.tsv", sample=SAMPLES)
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