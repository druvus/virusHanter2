# post_processing.smk

# Import custom functions
from scripts.functions import (
    read_file_as_blob,
    parse_bwa_flagstat,
    parse_fastp,
    plot_flagstat,
    plot_kaiju,
    plot_kraken,
    plot_blastn,
    alignment_stats,
    panel_report,
)

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
        "envs/bwa.yaml"
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
        shell("mkdir -p {os.path.dirname(index_prefix)}")
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
        coverage_plots_dir=f"{RESULT_FOLDER}/{{sample}}/COVERAGE_PLOTS",
    params:
        threshold=PLOT_THRESHOLD,
        num_refs=NUMBER_OF_PLOTS,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bam2plot.log"
    conda:
        "envs/bam2plot.yaml"
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

# Rule: Generate interactive report
rule generate_report:
    input:
        flagstat=rules.remove_host.output.flagstat,
        blastn_csv=rules.merge_checkv_blastn.output.merged_csv,
        kraken_csv=rules.wrangle_kraken.output.kraken_csv,
        kaiju_table=rules.kaiju_to_table.output.kaiju_table,
    output:
        report_html=f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html",
    params:
        result_folder=lambda wildcards: f"{RESULT_FOLDER}/{wildcards.sample}",
        secondary_host_name=SECONDARY_HOST_NAME,
    conda:
        "envs/panel.yaml"
    run:
        # Generate the report using the panel_report function
        report = panel_report(
            result_folder=params.result_folder,
            blastn_file=input.blastn_csv,
            kraken_file=input.kraken_csv,
            kaiju_table=input.kaiju_table,
            secondary_host=params.secondary_host_name,
        )
        # Save the report as an HTML file
        report.save(output.report_html, title=f"Report of {wildcards.sample}")

# Rule: Aggregate run information across samples
rule aggregate_run_information:
    input:
        reports=expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
    output:
        run_info_csv=f"{RESULT_FOLDER}/run_information_{Path(config['SAMPLES']).name}.csv",
    params:
        results_folder=RESULT_FOLDER,
    conda:
        "envs/panel.yaml"
    run:
        import pandas as pd
        from pathlib import Path

        def aggregate_sample_info(sample_folder):
            sample_folder = Path(sample_folder)
            sample_name = sample_folder.name

            # Read HTML report as blob
            report_html = read_file_as_blob(sample_folder / "REPORT" / f"{sample_name}.html")

            # Parse FASTP report
            fastp_report = next(sample_folder.rglob("FASTP/*.html"), None)
            fastp_df = parse_fastp(str(fastp_report)) if fastp_report else pd.DataFrame()

            # Extract relevant information
            sequencing_length = fastp_df.loc[fastp_df['description'] == 'Read1 Length', 'value'].values[0]
            number_reads = fastp_df.loc[fastp_df['description'] == 'Total Reads', 'value'].values[0]

            # Parse BWA flagstat
            flagstat_file = sample_folder / "logs" / "human_contamination_flagstat.txt"
            total_reads, percent_mapped = parse_bwa_flagstat(str(flagstat_file))

            # Read Kraken2 results
            kraken_csv = sample_folder / "KRAKEN" / f"{sample_name}.kraken.csv"
            kraken_df = pd.read_csv(kraken_csv)
            kraken_virus_percent = kraken_df.loc[kraken_df['domain'] == 'Viruses', 'percent'].sum()

            # Read Kaiju results
            kaiju_table = sample_folder / "KAIJU" / f"{sample_name}.kaiju.table.tsv"
            kaiju_df = pd.read_csv(kaiju_table, sep='\t')
            kaiju_virus_percent = kaiju_df['percent'].sum()
            top_virus_kaiju = '||'.join(kaiju_df['taxon_name'].head(10).tolist())

            # Read BLASTN results
            blastn_csv = sample_folder / "BLASTN" / f"{sample_name}.contigs.blastn.csv"
            blastn_df = pd.read_csv(blastn_csv)
            number_contigs = len(blastn_df)
            top_contigs_blastn = '||'.join(blastn_df['match_name'].head(5).tolist())

            # Compile sample information
            sample_info = {
                'sample_name': sample_name,
                'read_length': sequencing_length,
                'number_reads': number_reads,
                'mapped_to_human_percent': percent_mapped,
                'kraken_virus_percent': kraken_virus_percent,
                'kaiju_virus_percent': kaiju_virus_percent,
                'number_of_contigs': number_contigs,
                'top_contigs_blastn': top_contigs_blastn,
                'top_virus_kaiju': top_virus_kaiju,
                'report_html_blob': report_html,
            }

            return pd.DataFrame([sample_info])

        # Aggregate information for all samples
        samples_info = [aggregate_sample_info(Path(params.results_folder) / sample) for sample in SAMPLES]
        run_info_df = pd.concat(samples_info, ignore_index=True)

        # Save aggregated run information
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