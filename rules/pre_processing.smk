# pre_processing.smk

# Import custom functions
from scripts.functions import (
    paired_reads,
    common_suffix,
)

# Set variables
THREADS = config["THREADS"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(config["SAMPLES"]).name)
SAMPLES_FOLDER = config["SAMPLES"]
SUFFIX = common_suffix(SAMPLES_FOLDER)
SAMPLES = paired_reads(SAMPLES_FOLDER)
HUMAN_INDEX = config["HUMAN_INDEX"]
SECONDARY_HOST_INDEX = config.get("SECONDARY_HOST_INDEX", "")
SECONDARY_HOST_NAME = config.get("SECONDARY_HOST_NAME", "")
SECONDARY_HOST_OR_NOT = bool(SECONDARY_HOST_INDEX)

# Rule: Quality control with fastp
rule fastp:
    input:
        r1=lambda wildcards: f"{SAMPLES_FOLDER}/{wildcards.sample}1{SUFFIX}",
        r2=lambda wildcards: f"{SAMPLES_FOLDER}/{wildcards.sample}2{SUFFIX}",
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}_r1_trimmed.fq",
        r2=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}_r2_trimmed.fq",
        html_report=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}.fastp.html",
        json_report=f"{RESULT_FOLDER}/{{sample}}/FASTP/{{sample}}.fastp.json",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/fastp.log"
    threads: THREADS
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --report_title {wildcards.sample} \
            --thread {threads} \
            --html {output.html_report} \
            --json {output.json_report} \
            > {log} 2>&1
        """

# Rule: Align reads to human genome to identify contamination
rule bwa_human:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    params:
        index=HUMAN_INDEX,
    output:
        mapped_bam=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human.bam",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_human.log"
    threads: THREADS
    conda:
        "../envs/bwa.yaml"
    shell:
        # `-k 26` raises bwa-mem's minimum seed length from 19 to 26, matching
        # the original virusHanter setting for human contamination removal.
        # A higher seed length trades a little sensitivity for fewer spurious
        # short-seed matches against the human reference.
        """
        bwa mem -t {threads} -k 26 {params.index} {input.r1} {input.r2} | samtools sort -o {output.mapped_bam} - > {log} 2>&1
        """

# Rule: Remove reads mapped to human genome
rule remove_host:
    input:
        mapped_bam=rules.bwa_human.output.mapped_bam,
    output:
        unmapped_bam=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human_unmapped.bam",
        flagstat=f"{RESULT_FOLDER}/{{sample}}/logs/human_contamination_flagstat.txt",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/remove_host.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat {input.mapped_bam} > {output.flagstat}
        samtools view -b -f 12 {input.mapped_bam} > {output.unmapped_bam} 2>> {log}
        """

# Rule: Convert BAM to FASTQ after human host removal
rule bam_to_fastq_human:
    input:
        unmapped_bam=rules.remove_host.output.unmapped_bam,
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human_unmapped_r1.fastq",
        r2=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human_unmapped_r2.fastq",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bam_to_fastq_human.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} {input.unmapped_bam} > {log} 2>&1
        """

rule bwa_secondary_host:
    input:
        r1=rules.bam_to_fastq_human.output.r1,
        r2=rules.bam_to_fastq_human.output.r2,
    params:
        index=SECONDARY_HOST_INDEX,
        is_secondary_host=SECONDARY_HOST_OR_NOT,
    output:
        mapped_bam=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_secondary.bam",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bwa_secondary.log"
    threads: THREADS
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        if [ "{params.is_secondary_host}" = "True" ]; then
            bwa mem -t {threads} {params.index} {input.r1} {input.r2} 2> {log} | \
            samtools sort -o {output.mapped_bam} -
        else
            touch {output.mapped_bam}
        fi
        """

# Rule: Remove reads mapped to secondary host genome (optional)
rule remove_secondary_host:
    input:
        mapped_bam=rules.bwa_secondary_host.output.mapped_bam,
        unmapped_bam_human=rules.remove_host.output.unmapped_bam,
    output:
        unmapped_bam=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_secondary_unmapped.bam",
        flagstat=f"{RESULT_FOLDER}/{{sample}}/logs/secondary_contamination_flagstat.txt",
    params:
        is_secondary_host=str(SECONDARY_HOST_OR_NOT),
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/remove_secondary_host.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        if [ "{params.is_secondary_host}" = "True" ]; then
            samtools flagstat {input.mapped_bam} > {output.flagstat}
            samtools view -b -f 12 {input.mapped_bam} > {output.unmapped_bam} 2>> {log}
        else
            # Plain copy rather than a symlink so the output is portable
            # across machines and archivable in place.
            cp {input.unmapped_bam_human} {output.unmapped_bam}
            : > {output.flagstat}
        fi
        """

# Rule: Convert BAM to FASTQ after secondary host removal
rule bam_to_fastq_secondary:
    input:
        unmapped_bam=rules.remove_secondary_host.output.unmapped_bam,
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_final_unmapped_r1.fastq",
        r2=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_final_unmapped_r2.fastq",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/bam_to_fastq_secondary.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} {input.unmapped_bam} > {log} 2>&1
        """