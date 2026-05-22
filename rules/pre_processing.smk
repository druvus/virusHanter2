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
# Off by default to preserve byte-identical parity with the original
# virusHanter outputs. When TRUE, `remove_host` reads the markdup BAM
# and excludes flag-1024 reads from the host-removed FASTQs that feed
# MEGAHIT and the BWA-to-Kraken-hits coverage step.
DEDUPLICATE = config.get("DEDUPLICATE", "FALSE") == "TRUE"

# Host-removal backend. "bwa" is the parity default (bwa mem -k 26 +
# samtools); "hostile" calls Bede et al.'s hostile against the
# bundled T2T-CHM13 reference. Helpers below pick the output path
# from whichever backend is active; the lambda inputs scattered
# across the other .smk files use them rather than hard-coding
# `rules.bam_to_fastq_human.output`.
HOST_REMOVAL = config.get("HOST_REMOVAL", "bwa")
HOSTILE_INDEX = config.get("HOSTILE_INDEX", "")


def host_removed_r1(wildcards):
    """Resolve the FASTQ path for the host-removed R1 stream.

    Honours both `SECONDARY_HOST_OR_NOT` (which adds a second host
    layer) and `HOST_REMOVAL` (bwa vs hostile for the primary
    host). Used by every downstream consumer that previously
    referenced `rules.bam_to_fastq_human.output.r1` directly.
    """
    if SECONDARY_HOST_OR_NOT:
        return rules.bam_to_fastq_secondary.output.r1
    if HOST_REMOVAL == "hostile":
        return rules.hostile_human.output.r1
    return rules.bam_to_fastq_human.output.r1


def host_removed_r2(wildcards):
    if SECONDARY_HOST_OR_NOT:
        return rules.bam_to_fastq_secondary.output.r2
    if HOST_REMOVAL == "hostile":
        return rules.hostile_human.output.r2
    return rules.bam_to_fastq_human.output.r2


def host_flagstat(wildcards=None):
    """Resolve the path of the primary-host flagstat artefact.

    Both backends emit a flagstat-format text file at the same
    canonical path under `logs/`; the bwa backend gets it from
    `samtools flagstat` directly, the hostile backend uses hostile's
    own stats output reformatted to the flagstat shape.
    """
    if HOST_REMOVAL == "hostile":
        return rules.hostile_human.output.flagstat
    return rules.remove_host.output.flagstat

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

# Rule: Mark PCR duplicates on the host-aligned BAM.
#
# The textual stats file is always emitted (consumed by MultiQC and
# reporthanter). The markdup BAM is emitted as well so `remove_host`
# can read it when `DEDUPLICATE: TRUE` is set in config. With
# `DEDUPLICATE: FALSE` (the default), `remove_host` ignores this BAM
# and continues to consume the un-marked bwa_human output — the
# pipeline is then byte-identical to the pre-markdup behaviour.
#
# `samtools markdup` requires MC/MS tags added by `fixmate`, which in
# turn needs a name-sorted BAM. `bwa_human` writes a coordinate-sorted
# BAM, so the chain is: name-sort -> fixmate -> coord-sort -> markdup.
rule markdup_human:
    input:
        mapped_bam=rules.bwa_human.output.mapped_bam,
    output:
        markdup_bam=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human_markdup.bam",
        stats=f"{RESULT_FOLDER}/{{sample}}/logs/human_markdup_stats.txt",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/markdup_human.log"
    threads: THREADS
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -n -@ {threads} -O bam {input.mapped_bam} 2>> {log} \
          | samtools fixmate -m -@ {threads} -O bam - - 2>> {log} \
          | samtools sort -@ {threads} -O bam - 2>> {log} \
          | samtools markdup -@ {threads} -s -f {output.stats} - {output.markdup_bam} 2>> {log}
        """


# Rule: Remove reads mapped to human genome.
#
# When `DEDUPLICATE: TRUE`, this rule reads the markdup BAM and
# excludes flag-1024 reads so downstream assembly and coverage are
# not inflated by PCR duplicates. The default (`DEDUPLICATE: FALSE`)
# reads the un-marked bwa_human BAM, preserving byte-identical parity
# with the original virusHanter pipeline.
rule remove_host:
    input:
        mapped_bam=(
            rules.markdup_human.output.markdup_bam
            if DEDUPLICATE
            else rules.bwa_human.output.mapped_bam
        ),
    output:
        unmapped_bam=f"{RESULT_FOLDER}/{{sample}}/bwa/{{sample}}_human_unmapped.bam",
        flagstat=f"{RESULT_FOLDER}/{{sample}}/logs/human_contamination_flagstat.txt",
    params:
        dup_filter="-F 1024" if DEDUPLICATE else "",
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/remove_host.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat {input.mapped_bam} > {output.flagstat}
        samtools view -b -f 12 {params.dup_filter} {input.mapped_bam} > {output.unmapped_bam} 2>> {log}
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

# Rule: Host removal via hostile.
#
# Alternative to the bwa_human / markdup_human / remove_host /
# bam_to_fastq_human chain when `HOST_REMOVAL: "hostile"` is set in
# config. hostile (Bede et al.) wraps minimap2 against a bundled
# T2T-CHM13 + alt + decoy human reference; it catches host reads
# from telomeres, centromeres and segmental duplications that
# GRCh38-based BWA misses. Outputs the same R1 / R2 paired FASTQ
# pair that bam_to_fastq_human produces (different path so both
# backends can coexist on disk) plus a `flagstat`-shape stats file
# the report layer already knows how to parse.
#
# `HOSTILE_INDEX` is the directory hostile downloads the T2T-CHM13
# bundle to (see hostile docs). Leave empty to let hostile manage
# its own cache.
rule hostile_human:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    output:
        r1=f"{RESULT_FOLDER}/{{sample}}/hostile/{{sample}}_human_unmapped_r1.fastq",
        r2=f"{RESULT_FOLDER}/{{sample}}/hostile/{{sample}}_human_unmapped_r2.fastq",
        # Distinct path from the BWA chain's flagstat so both rule
        # families can co-exist in the workflow at parse time.
        # `host_flagstat()` picks the right path based on
        # HOST_REMOVAL.
        flagstat=f"{RESULT_FOLDER}/{{sample}}/logs/hostile_contamination_flagstat.txt",
        stats_json=f"{RESULT_FOLDER}/{{sample}}/hostile/hostile_stats.json",
    params:
        out_dir=f"{RESULT_FOLDER}/{{sample}}/hostile",
        # Default to the viral-mask + phage-mask T2T-CHM13 index.
        # The bare `human-t2t-hla` index would filter reads from
        # endogenous retroviruses, phages and other host-embedded
        # viral elements as host — exactly the reads a viral
        # metagenomics pipeline must preserve. The masked variant
        # leaves the human genome intact except for those regions
        # that align to a curated viral / phage reference, so the
        # filter remains specific to the true host backbone.
        index_arg=lambda wildcards: (
            f"--index {HOSTILE_INDEX}"
            if HOSTILE_INDEX
            else "--index human-t2t-hla.rs-viral-202401_ml-phage-202401"
        ),
    threads: THREADS
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/hostile_human.log"
    conda:
        "../envs/hostile.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        # hostile clean: write the host-removed FASTQs to {params.out_dir}
        # and emit a JSON stats sidecar. Use the short-read minimap2
        # preset (the default for paired Illumina input).
        hostile clean \
            --fastq1 {input.r1} --fastq2 {input.r2} \
            --output {params.out_dir} \
            --threads {threads} \
            {params.index_arg} \
            > {output.stats_json} 2> {log}

        # hostile names its outputs after the input filenames. The
        # exact basenames depend on the input; rename / move so the
        # rule's declared outputs always sit at the expected paths.
        mv {params.out_dir}/*.clean_1.fastq.gz {output.r1}.gz 2>> {log} || \
            mv {params.out_dir}/*.clean_1.fastq {output.r1} 2>> {log}
        mv {params.out_dir}/*.clean_2.fastq.gz {output.r2}.gz 2>> {log} || \
            mv {params.out_dir}/*.clean_2.fastq {output.r2} 2>> {log}
        # Decompress if the output landed gzipped.
        if [ -f {output.r1}.gz ]; then gunzip {output.r1}.gz; fi
        if [ -f {output.r2}.gz ]; then gunzip {output.r2}.gz; fi

        # Reformat hostile's JSON into samtools-flagstat shape so
        # the downstream consumers (FlagstatProcessor in reporthanter
        # and aggregate_run_information) keep working unchanged.
        # hostile's JSON reports `reads_in` and `reads_out`; flagstat
        # convention is "paired in sequencing" and "with itself and
        # mate mapped" referring to the host-aligned pairs. The
        # newline characters are doubly escaped so the outer
        # Snakefile (Python) parse leaves them as `\\n`, bash passes
        # them through verbatim, and the inner `python -c "..."` parse
        # finally interprets them as real newlines.
        python -c "
import json
stats = json.load(open('{output.stats_json}'))
s = stats[0] if isinstance(stats, list) else stats
reads_in = int(s.get('reads_in', 0))
reads_out = int(s.get('reads_out', 0))
mapped = reads_in - reads_out
with open('{output.flagstat}', 'w') as fh:
    fh.write(str(reads_in) + ' + 0 paired in sequencing\\n')
    fh.write(str(mapped) + ' + 0 with itself and mate mapped\\n')
" >> {log} 2>&1
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