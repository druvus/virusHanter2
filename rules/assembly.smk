# assembly.smk
#
# Assembly, polishing, contig annotation, and contamination assessment.
# Ported from the original virusHanter/Snakefile (rules: megahit, pilon,
# wrangle_pilon, blastn, checkv, merge_checkv_blastn).

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
        "envs/megahit.yaml"
    run:
        # MEGAHIT refuses to run if the output directory already exists.
        shell("rm -rf {params.out_dir}")
        shell(
            "megahit "
            "-1 {input.r1} -2 {input.r2} "
            "-o {params.out_dir} "
            "--out-prefix {wildcards.sample} "
            "-t {threads} "
            "2> {log}"
        )
        # Drop intermediate MEGAHIT files; keep the contigs FASTA.
        shell("ls -d -1 {params.out_dir}/* | grep -v .fa | xargs rm -rf")

        # If MEGAHIT produced no contigs, emit a dummy contig so downstream
        # rules (BLASTN, CheckV, Pilon) still have an input. Mirrors the
        # original virusHanter behavior.
        if Path(output.contigs).read_text() == "":
            with open(output.contigs, "w") as f:
                f.write(">DUMMY_CONTIG\n")
                f.write("TTAACCTTGG" * 20 + "\n")


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
        # mem_mb is derived from PILON_MEM so the scheduler and the JVM Xmx
        # flag stay in lockstep when PILON_MEM is tuned.
        mem_mb=lambda wildcards, params: int(params.pilon_mem.rstrip("Gg")) * 1024,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/pilon.log",
    conda:
        "envs/pilon.yaml"
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
        "envs/panel.yaml"
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
        "envs/blastn.yaml"
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
        "envs/checkv.yaml"
    shell:
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
        "envs/panel.yaml"
    run:
        import pandas as pd

        blastn_df = pd.read_csv(input.blastn)
        checkv_df = (
            pd.read_csv(input.checkv, sep="\t")
            .rename(columns={"contig_id": "name"})
            [["name", "total_genes", "viral_genes", "host_genes", "provirus"]]
        )
        merged = pd.merge(blastn_df, checkv_df, on="name", how="left")
        merged.to_csv(output.merged_csv, index=False)
