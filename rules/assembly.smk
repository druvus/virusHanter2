# assembly.smk
#
# Assembly, polishing, contig annotation, and contamination assessment.
# Ported from the original virusHanter/Snakefile (rules: megahit, pilon,
# wrangle_pilon, blastn, checkv, merge_checkv_blastn).

import platform

from scripts.functions import (
    fastx_file_to_df,
    run_blastn,
)

# Set variables (RESULT_FOLDER and THREADS are also defined in the other
# included rule files; Snakemake tolerates redefinition with identical values).
THREADS = config["THREADS"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(config["SAMPLES"]).name)
CONTIG_LENGTH = config.get("CONTIG_LENGTH", 500)
BLASTN_DB = config["BLASTN_DB"]
CHECKV_DB = config["CHECKV_DB"]
PILON_MEM = config.get("PILON_MEM", "50G")
PILON_MEM_MB = int(PILON_MEM.rstrip("Gg")) * 1024

# Optional geNomad classifier (off by default; turn on with
# `GENOMAD: "TRUE"` in config plus a populated `GENOMAD_DB`).
RUN_GENOMAD = config.get("GENOMAD", "FALSE") == "TRUE"
GENOMAD_DB = config.get("GENOMAD_DB", "")

# Optional QUAST assembly assessment (off by default). When TRUE,
# `rule quast_megahit` runs against each sample's MEGAHIT contigs and
# its report dir is fed to MultiQC for batch-level QC.
RUN_QUAST = config.get("QUAST", "FALSE") == "TRUE"


# Rule: De novo assembly with MEGAHIT
rule megahit:
    input:
        r1=lambda wildcards: rules.bam_to_fastq_human.output.r1 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r1,
        r2=lambda wildcards: rules.bam_to_fastq_human.output.r2 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r2,
    output:
        contigs=f"{RESULT_FOLDER}/{{sample}}/MEGAHIT/{{sample}}.contigs.fa",
    params:
        out_dir=f"{RESULT_FOLDER}/{{sample}}/MEGAHIT",
    threads: THREADS
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/megahit.log",
    conda:
        "../envs/megahit.yaml"
    run:
        import subprocess

        # MEGAHIT refuses to run if the output directory already exists.
        shell("rm -rf {params.out_dir}")
        # MEGAHIT defaults to allocating 90% of detected RAM up front,
        # which SIGSEGVs on a memory-tight host (e.g. an 18 GB laptop with
        # other processes already running). Bound the request to a
        # configurable fraction so it fits on small hosts; on a Linux box
        # with abundant RAM this is still a sensible cap.
        mem_fraction = float(config.get("MEGAHIT_MEM_FRACTION", 0.5))
        # Apple Silicon bioconda osx-arm64 megahit has three reproducible
        # quirks: (1) `megahit_core_popcnt` segfaults on `count -k 21`
        # regardless of memory; (2) `megahit_core_no_hw_accel` segfaults
        # at `-t > 2`; (3) `megahit_core_no_hw_accel count -k 21` also
        # SIGSEGVs on small inputs (the smoke fixture, e.g.) — raising
        # the minimum k past 21 avoids the buggy path. Force the
        # no_hw_accel variant, cap threads at 2, and start at k=27 on
        # Darwin/arm64. All three flags are no-ops or harmless on Linux,
        # but k=27 is a small assembly-quality concession so it is gated
        # on the platform check rather than applied globally.
        is_apple_silicon = platform.system() == "Darwin" and platform.machine() == "arm64"
        no_hw_accel = "--no-hw-accel " if is_apple_silicon else ""
        kmin_flag = "--k-min 27 " if is_apple_silicon else ""
        mh_threads = min(threads, 2) if is_apple_silicon else threads
        try:
            shell(
                "megahit "
                "-1 {input.r1} -2 {input.r2} "
                "-o {params.out_dir} "
                "--out-prefix {wildcards.sample} "
                f"-t {mh_threads} "
                f"-m {mem_fraction} "
                f"{no_hw_accel}"
                f"{kmin_flag}"
                "2> {log}"
            )
        except subprocess.CalledProcessError:
            # MEGAHIT can SIGSEGV on tiny inputs. Make sure the output file
            # exists so the dummy-contig fallback below writes to it.
            Path(params.out_dir).mkdir(parents=True, exist_ok=True)
            Path(output.contigs).touch()

        # Drop intermediate MEGAHIT files; keep the contigs FASTA. The grep
        # may exit non-zero when the directory only contains the .fa output;
        # tolerate that.
        shell(
            "ls -d -1 {params.out_dir}/* 2>/dev/null "
            "| grep -v .fa | xargs rm -rf || true"
        )

        # If MEGAHIT produced no contigs (or crashed), emit a dummy contig
        # so downstream rules (BLASTN, CheckV, Pilon) still have an input.
        # Mirrors the original virusHanter behavior.
        if Path(output.contigs).read_text() == "":
            with open(output.contigs, "w") as f:
                f.write(">DUMMY_CONTIG\n")
                f.write("TTAACCTTGG" * 20 + "\n")


# Rule: QUAST assembly assessment on the raw MEGAHIT contigs.
#
# Reports N50, largest contig, total assembled length, GC% and other
# standard assembly metrics. Runs on the un-polished MEGAHIT output so
# the metrics describe the assembler's behaviour directly; Pilon
# improvements are a separate concern. Output is consumed by MultiQC.
rule quast_megahit:
    input:
        contigs=rules.megahit.output.contigs,
    output:
        report_dir=directory(f"{RESULT_FOLDER}/{{sample}}/QUAST"),
        report_tsv=f"{RESULT_FOLDER}/{{sample}}/QUAST/report.tsv",
    threads: 2
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/quast.log",
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast.py \
            --threads {threads} \
            --output-dir {output.report_dir} \
            --labels {wildcards.sample} \
            {input.contigs} \
            > {log} 2>&1
        """


# Rule: Polish contigs with Pilon
rule pilon:
    input:
        contigs=rules.megahit.output.contigs,
        r1=lambda wildcards: rules.bam_to_fastq_human.output.r1 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r1,
        r2=lambda wildcards: rules.bam_to_fastq_human.output.r2 if not SECONDARY_HOST_OR_NOT else rules.bam_to_fastq_secondary.output.r2,
    output:
        contigs_bam=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_contigs.bam",
        improved_contigs=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}_improved_contigs.fasta",
    params:
        index_folder=f"{RESULT_FOLDER}/{{sample}}/PILON/bwa",
        pilon_folder=f"{RESULT_FOLDER}/{{sample}}/PILON",
        pilon_mem=PILON_MEM,
    threads: THREADS
    resources:
        # PILON_MEM_MB is derived from PILON_MEM at workflow-parse time so the
        # scheduler and the JVM Xmx flag stay in lockstep.
        mem_mb=PILON_MEM_MB,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/pilon.log",
    conda:
        "../envs/pilon.yaml"
    run:
        # Rebuild a temporary BWA index of the assembly.
        shell("rm -rf {params.index_folder}")
        shell("mkdir -p {params.index_folder}")
        index_prefix = f"{params.index_folder}/{wildcards.sample}"
        shell("bwa index -p {index_prefix} {input.contigs} >> {log} 2>&1")

        # Map reads back to contigs and sort.
        shell(
            "bwa mem -t {threads} {index_prefix} {input.r1} {input.r2} 2>> {log} "
            "| samtools view -h -O bam "
            "| samtools sort -o {output.contigs_bam}"
        )
        shell("samtools index {output.contigs_bam}")

        # Polish.
        shell(
            "pilon -Xmx{params.pilon_mem} --threads {threads} "
            "--genome {input.contigs} --frags {output.contigs_bam} "
            "--outdir {params.pilon_folder} "
            "--output {wildcards.sample}_improved_contigs "
            ">> {log} 2>&1"
        )

        # Drop the temporary index.
        shell("rm -rf {params.index_folder}")


# Rule: Convert polished contigs FASTA into a length-filtered CSV
rule wrangle_pilon:
    input:
        contigs=rules.pilon.output.improved_contigs,
    output:
        csv=f"{RESULT_FOLDER}/{{sample}}/PILON/{{sample}}.contigs.csv",
    params:
        min_len=CONTIG_LENGTH,
    conda:
        "../envs/panel.yaml"
    run:
        df = fastx_file_to_df(input.contigs)
        df = df.assign(sample_id=wildcards.sample)
        df = df.loc[lambda x: x.read_len > params.min_len]
        df.to_csv(output.csv, index=False)


# Rule: Annotate contigs with BLASTN
rule blastn:
    input:
        contigs=rules.wrangle_pilon.output.csv,
    output:
        blast=f"{RESULT_FOLDER}/{{sample}}/BLASTN/{{sample}}.contigs.blastn.csv",
    params:
        temp_file=f"{RESULT_FOLDER}/{{sample}}/BLASTN/temp.blastn.fasta",
        db=BLASTN_DB,
    threads: THREADS
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/blastn.log",
    conda:
        "../envs/blastn.yaml"
    run:
        Path(output.blast).parent.mkdir(parents=True, exist_ok=True)
        df = run_blastn(
            contigs_csv=input.contigs,
            db=params.db,
            temp_file=params.temp_file,
            threads=threads,
        )
        df.to_csv(output.blast, index=False)

        if Path(params.temp_file).exists():
            os.remove(params.temp_file)


# Rule: Assess host contamination of contigs with CheckV
rule checkv:
    input:
        contigs=rules.pilon.output.improved_contigs,
    output:
        checkv=f"{RESULT_FOLDER}/{{sample}}/CHECKV/{{sample}}.contamination.tsv",
    params:
        db=CHECKV_DB,
        folder=f"{RESULT_FOLDER}/{{sample}}/CHECKV",
    threads: THREADS
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/checkv.log",
    conda:
        "../envs/checkv.yaml"
    shell:
        # Known issue on macOS (both osx-64 via Rosetta and native osx-arm64):
        # CheckV 1.0.3's `search_hmms` reports "80 hmmsearch tasks failed"
        # even when every .hmmout file ends with the `# [ok]` marker, because
        # sp.Popen(cmd, shell=True).wait() returns non-zero in that
        # multiprocessing.Pool worker context. The bug is independent of the
        # number of threads (also fails at -t 1) and the hmmer build, and
        # there is no newer CheckV release on bioconda. Production / Phase 6
        # parity runs must use Linux.
        """
        checkv contamination \
            -d {params.db} \
            {input.contigs} \
            {params.folder} \
            -t {threads} \
            2> {log}

        mv {params.folder}/contamination.tsv {output.checkv}
        # Drop CheckV intermediates; keep the contamination TSV only.
        ls -d -1 {params.folder}/* | grep -v .tsv | xargs rm -rf
        """


# Rule: Merge CheckV contamination calls into the BLASTN annotation table
rule merge_checkv_blastn:
    input:
        checkv=rules.checkv.output.checkv,
        blastn=rules.blastn.output.blast,
    output:
        merged_csv=f"{RESULT_FOLDER}/{{sample}}/CHECKV/{{sample}}.merged.csv",
    conda:
        "../envs/panel.yaml"
    run:
        import pandas as pd

        blastn_df = pd.read_csv(input.blastn)
        checkv_df = (
            pd.read_csv(input.checkv, sep="\t")
            .rename(columns={"contig_id": "name"})
            [["name", "total_genes", "viral_genes", "host_genes", "provirus"]]
        )
        # Inner join on `name` matches the original virusHanter behaviour:
        # contigs without a CheckV entry are dropped from the merged table.
        merged = pd.merge(blastn_df, checkv_df, on="name", how="inner")
        merged.to_csv(output.merged_csv, index=False)


# Rule: Optional second viral-contig classifier (geNomad).
#
# Off by default. Turn on with `GENOMAD: "TRUE"` and a populated
# `GENOMAD_DB` in config. Stores a summary TSV under GENOMAD/ alongside
# the existing CheckV outputs. Does NOT feed the per-sample HTML or
# the run-aggregation CSV; the merged_csv schema and every other
# parity-locked output stay byte-identical.
#
# Pulls the Pilon-polished contigs (same input CheckV uses) and runs
# `genomad end-to-end` to classify each contig as plasmid, viral, or
# chromosomal. The summary file is the per-sample headline output:
# `<sample>_summary/<sample>_virus_summary.tsv`.
rule genomad:
    input:
        contigs=rules.pilon.output.improved_contigs,
    output:
        summary=f"{RESULT_FOLDER}/{{sample}}/GENOMAD/{{sample}}_summary/{{sample}}_virus_summary.tsv",
    params:
        db=GENOMAD_DB,
        out_dir=f"{RESULT_FOLDER}/{{sample}}/GENOMAD",
    threads: THREADS
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/genomad.log"
    conda:
        "../envs/genomad.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        genomad end-to-end \
            --threads {threads} \
            --cleanup \
            {input.contigs} \
            {params.out_dir} \
            {params.db} \
            > {log} 2>&1
        """
