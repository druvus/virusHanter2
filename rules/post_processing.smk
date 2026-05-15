# post_processing.smk
#
# Coverage alignment, report generation (via the reporthanter CLI), run
# aggregation, and optional cleanup.

from scripts.functions import read_file_as_blob

# Set variables
THREADS = config["THREADS"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(config["SAMPLES"]).name)
VIRUS_PARQUET = config["VIRUS_PARQUET"]
PLOT_THRESHOLD = config["PLOT_THRESHOLD"]
NUMBER_OF_PLOTS = config["NUMBER_OF_PLOTS"]
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

        # Read top Kraken2 viral taxa
        kraken_df = pd.read_csv(input.kraken_csv)
        top_tax_ids = kraken_df.loc[kraken_df.domain == "Viruses", "taxonomy_id"].head(20).tolist()

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

# Rule: Generate coverage plots
rule bam2plot:
    input:
        bam=rules.bwa_align_to_kraken_hits.output.bam,
    output:
        coverage_plots_dir=directory(f"{RESULT_FOLDER}/{{sample}}/COVERAGE_PLOTS"),
    params:
        threshold=PLOT_THRESHOLD,
        num_refs=NUMBER_OF_PLOTS,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bam2plot.log"
    conda:
        "../envs/bam2plot.yaml"
    shell:
        """
        mkdir -p {output.coverage_plots_dir}
        bam2plot from_bam \
            -b {input.bam} \
            -o {output.coverage_plots_dir} \
            -t {params.threshold} \
            -p svg \
            -n {params.num_refs} \
            > {log} 2>&1
        """

# Rule: Generate interactive report via the reporthanter CLI
#
# The reporthanter package is expected to be installed in the driver conda env
# (virushanter) that runs Snakemake, so this rule deliberately omits a
# `conda:` directive.
rule generate_report:
    input:
        flagstat=rules.remove_host.output.flagstat,
        secondary_flagstat=rules.remove_secondary_host.output.flagstat,
        fastp_json=rules.fastp.output.json_report,
        blastn_csv=rules.merge_checkv_blastn.output.merged_csv,
        kraken_csv=rules.wrangle_kraken.output.kraken_csv,
        kaiju_table=rules.kaiju_to_table.output.kaiju_table,
        coverage_dir=rules.bam2plot.output.coverage_plots_dir,
    output:
        report_html=f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html",
    params:
        secondary_args=(
            lambda wildcards: (
                f"--secondary_flagstat_file {RESULT_FOLDER}/{wildcards.sample}/logs/secondary_contamination_flagstat.txt "
                f"--secondary_host {SECONDARY_HOST_NAME}"
            )
            if SECONDARY_HOST_OR_NOT
            else ""
        ),
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/reporthanter.log",
    shell:
        """
        reporthanter \
            --blastn_file {input.blastn_csv} \
            --kraken_file {input.kraken_csv} \
            --kaiju_table {input.kaiju_table} \
            --fastp_json {input.fastp_json} \
            --flagstat_file {input.flagstat} \
            --coverage_folder {input.coverage_dir} \
            --output {output.report_html} \
            --sample_name {wildcards.sample} \
            {params.secondary_args} \
            > {log} 2>&1
        """

# Rule: Aggregate run information across samples
#
# Runs in the driver conda env (virushanter) which has reporthanter + pandas
# installed, so no per-rule `conda:` directive is needed. Parsing of fastp
# JSON and BWA flagstat output is delegated to the reportHanter processors to
# keep a single implementation across the pipeline and the HTML report.
rule aggregate_run_information:
    input:
        reports=expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
    output:
        run_info_csv=f"{RESULT_FOLDER}/run_information_{Path(config['SAMPLES']).name}.csv",
    params:
        results_folder=RESULT_FOLDER,
    log:
        f"{RESULT_FOLDER}/logs/aggregate_run_information.log",
    run:
        import json
        import pandas as pd
        from pathlib import Path

        from reporthanter import FlagstatProcessor

        flagstat_proc = FlagstatProcessor()

        def aggregate_sample_info(sample_folder: Path) -> pd.DataFrame:
            sample_name = sample_folder.name

            # Per-sample HTML, encoded as hex (matches the original virusHanter
            # run_information_<batch>.csv schema).
            report_html = read_file_as_blob(
                sample_folder / "REPORT" / f"{sample_name}.html"
            )

            # fastp summary (JSON) — read directly; only a couple of fields
            # are needed here, so a full FastpProcessor pass is unnecessary.
            fastp_json_path = sample_folder / "FASTP" / f"{sample_name}.fastp.json"
            with open(fastp_json_path) as fh:
                fastp_summary = json.load(fh).get("summary", {})
            before = fastp_summary.get("before_filtering", {})
            read_length = before.get("read1_mean_length", "")
            number_reads = before.get("total_reads", 0)

            # BWA flagstat against the human host.
            flagstat_path = sample_folder / "logs" / "human_contamination_flagstat.txt"
            flagstat_df = flagstat_proc.process(str(flagstat_path))
            flagstat_lookup = dict(zip(flagstat_df["metric"], flagstat_df["value"]))
            percent_mapped = flagstat_lookup.get("percent_mapped", 0.0)

            # Kraken2 / Kaiju percent-viral summaries.
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

            # BLASTN summary.
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

        samples_info = [
            aggregate_sample_info(Path(params.results_folder) / sample)
            for sample in SAMPLES
        ]
        run_info_df = pd.concat(samples_info, ignore_index=True)
        run_info_df.to_csv(output.run_info_csv, index=False)

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