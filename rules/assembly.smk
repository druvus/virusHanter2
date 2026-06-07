# assembly.smk
#
# Assembly, polishing, contig annotation, and contamination assessment.
# Ported from the original virusHanter/Snakefile (rules: megahit, pilon,
# wrangle_pilon, blastn, checkv, merge_checkv_blastn).
#
# Multi-assembler structure: every contig-producing rule below MEGAHIT
# and metaSPAdes carries an `{assembler}` wildcard so each downstream
# step (Pilon, BLASTN, CheckV, geNomad, QUAST) runs once per
# (sample, assembler) pair. The active assembler list lives in
# config[ASSEMBLERS] and is loaded in Snakefile as `ASSEMBLERS`.

import platform
from pathlib import Path

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

# Validate PILON_MEM: must be a positive integer followed by 'G' or 'g',
# and must not exceed available system RAM by an obvious margin.
import re as _re
import os as _os

_pilon_mem_match = _re.fullmatch(r"(\d+)[Gg]", str(PILON_MEM))
if not _pilon_mem_match:
    raise WorkflowError(
        f"config[PILON_MEM] must be a positive integer followed by 'G' "
        f"(e.g. '50G'), got {PILON_MEM!r}."
    )
_pilon_mem_gb = int(_pilon_mem_match.group(1))
if _pilon_mem_gb < 1:
    raise WorkflowError(
        f"config[PILON_MEM] must be at least '1G', got {PILON_MEM!r}."
    )
try:
    _available_ram_gb = int(
        _os.sysconf("SC_PAGE_SIZE") * _os.sysconf("SC_PHYS_PAGES") / (1024 ** 3)
    )
    if _pilon_mem_gb > _available_ram_gb:
        import sys as _sys
        print(
            f"WARNING: config[PILON_MEM] ({PILON_MEM}) exceeds available "
            f"system RAM ({_available_ram_gb} GB). Pilon may be killed by the "
            "OS OOM reaper. Reduce PILON_MEM or run on a larger host.",
            file=_sys.stderr,
        )
except (AttributeError, ValueError):
    pass  # sysconf not available on this platform; skip the RAM check.

PILON_MEM_MB = _pilon_mem_gb * 1024

# Optional geNomad classifier (off by default; turn on with
# `GENOMAD: "TRUE"` in config plus a populated `GENOMAD_DB`).
RUN_GENOMAD = config.get("GENOMAD", "FALSE") == "TRUE"
GENOMAD_DB = config.get("GENOMAD_DB", "")

# Validate GENOMAD_SPLITS: must be a non-negative integer (0 = mmseqs
# auto-split; any positive value partitions the mmseqs search to cap
# peak RAM). Validated here so a bad value raises at parse time rather
# than inside a running geNomad job.
_genomad_splits_raw = config.get("GENOMAD_SPLITS", 4)
try:
    _GENOMAD_SPLITS = int(_genomad_splits_raw)
except (TypeError, ValueError):
    raise WorkflowError(
        f"config[GENOMAD_SPLITS] must be a non-negative integer (0 = "
        f"mmseqs auto-split), got {_genomad_splits_raw!r}."
    )
if _GENOMAD_SPLITS < 0:
    raise WorkflowError(
        f"config[GENOMAD_SPLITS] must be >= 0, got {_GENOMAD_SPLITS}."
    )

# Optional QUAST assembly assessment (off by default). When TRUE,
# `rule quast_per_assembler` runs against each (sample, assembler) and
# its report dir is fed to MultiQC for batch-level QC.
RUN_QUAST = config.get("QUAST", "FALSE") == "TRUE"


def assembler_contigs(wildcards):
    """Map an {assembler} wildcard to its raw contigs FASTA.

    Lets downstream rules consume the contigs without caring which
    assembler produced them. Both megahit and metaspades land their
    output at the same path shape:
    `{sample}/{assembler}/{sample}.contigs.fa`.
    """
    return (
        f"{RESULT_FOLDER}/{wildcards.sample}/{wildcards.assembler}"
        f"/{wildcards.sample}.contigs.fa"
    )


# Rule: De novo assembly with MEGAHIT
rule megahit:
    input:
        r1=host_removed_r1,
        r2=host_removed_r2,
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
    script:
        "../scripts/run_megahit.py"


# Rule: De novo assembly with metaSPAdes
#
# Runs SPAdes in `--meta` mode on the same host-removed read pool that
# feeds MEGAHIT. metaSPAdes is markedly less crash-prone than MEGAHIT
# on Apple Silicon but can still exit non-zero on libraries it
# considers too small (it imposes a per-library minimum that MEGAHIT
# does not). The rule mirrors MEGAHIT's "dummy contig on failure"
# fallback so the per-assembler DAG stays uniform; downstream rules do
# not need to special-case an absent SPAdes output.
rule metaspades:
    input:
        r1=host_removed_r1,
        r2=host_removed_r2,
    output:
        contigs=f"{RESULT_FOLDER}/{{sample}}/metaSPAdes/{{sample}}.contigs.fa",
    params:
        out_dir=f"{RESULT_FOLDER}/{{sample}}/metaSPAdes",
        mode="meta",
    threads: THREADS
    resources:
        mem_mb=32000,
        runtime=360,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/metaspades.log",
    conda:
        "../envs/spades.yaml"
    script:
        "../scripts/run_spades.py"


# Rule: De novo assembly with rnaviralSPAdes
#
# SPAdes variant tuned for RNA viral libraries: handles the
# coverage variance and large insert-size distributions typical
# of host-depleted RNA virus samples better than `--meta` does.
# Mirrors the `metaspades` rule's command shape and dummy-contig
# fallback so the downstream {assembler}-wildcard chain stays
# uniform.
rule rnaviralspades:
    input:
        r1=host_removed_r1,
        r2=host_removed_r2,
    output:
        contigs=f"{RESULT_FOLDER}/{{sample}}/rnaviralSPAdes/{{sample}}.contigs.fa",
    params:
        out_dir=f"{RESULT_FOLDER}/{{sample}}/rnaviralSPAdes",
        mode="rnaviral",
    threads: THREADS
    resources:
        mem_mb=32000,
        runtime=360,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/rnaviralspades.log",
    conda:
        "../envs/spades.yaml"
    script:
        "../scripts/run_spades.py"


# Rule: QUAST assembly assessment on the raw assembler contigs.
#
# Reports N50, largest contig, total assembled length, GC% and other
# standard assembly metrics. Runs on the un-polished assembler output
# so the metrics describe the assembler's behaviour directly; Pilon
# improvements are a separate concern. One QUAST report per
# (sample, assembler); the report dir is consumed by MultiQC.
rule quast_per_assembler:
    input:
        contigs=assembler_contigs,
    output:
        report_dir=directory(f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/QUAST"),
        report_tsv=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/QUAST/report.tsv",
    threads: 2
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/quast_{{assembler}}.log",
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast.py \
            --threads {threads} \
            --output-dir {output.report_dir} \
            --labels {wildcards.sample}_{wildcards.assembler} \
            {input.contigs} \
            > {log} 2>&1
        """


# Rule: Polish contigs with Pilon (per assembler)
rule pilon:
    input:
        contigs=assembler_contigs,
        r1=host_removed_r1,
        r2=host_removed_r2,
    output:
        contigs_bam=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/PILON/{{sample}}_contigs.bam",
        improved_contigs=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/PILON/{{sample}}_improved_contigs.fasta",
    params:
        index_folder=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/PILON/bwa",
        pilon_folder=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/PILON",
        pilon_mem=PILON_MEM,
    threads: THREADS
    resources:
        # PILON_MEM_MB is derived from PILON_MEM at workflow-parse time so the
        # scheduler and the JVM Xmx flag stay in lockstep.
        mem_mb=PILON_MEM_MB,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/pilon_{{assembler}}.log",
    conda:
        "../envs/pilon.yaml"
    shell:
        """
        rm -rf {params.index_folder}
        mkdir -p {params.index_folder}
        INDEX_PREFIX={params.index_folder}/{wildcards.sample}
        bwa index -p $INDEX_PREFIX {input.contigs} >> {log} 2>&1
        bwa mem -t {threads} $INDEX_PREFIX {input.r1} {input.r2} 2>> {log} \
            | samtools view -h -O bam \
            | samtools sort -o {output.contigs_bam}
        samtools index {output.contigs_bam}
        pilon -Xmx{params.pilon_mem} --threads {threads} \
            --genome {input.contigs} --frags {output.contigs_bam} \
            --outdir {params.pilon_folder} \
            --output {wildcards.sample}_improved_contigs \
            >> {log} 2>&1
        rm -rf {params.index_folder}
        """


# Rule: Convert polished contigs FASTA into a length-filtered CSV
#
# Carries the {assembler} wildcard into a column on the CSV so every
# downstream consumer (BLASTN merge, per_virus_metrics, the report)
# knows which assembler produced each contig without re-deriving it
# from the file path.
rule wrangle_pilon:
    input:
        contigs=rules.pilon.output.improved_contigs,
    output:
        csv=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/PILON/{{sample}}.contigs.csv",
    params:
        min_len=CONTIG_LENGTH,
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/panel.yaml"
    run:
        df = fastx_file_to_df(input.contigs)
        df = df.assign(
            sample_id=wildcards.sample,
            assembler=wildcards.assembler,
        )
        df = df.loc[lambda x: x.read_len > params.min_len]
        df.to_csv(output.csv, index=False)


# Rule: Annotate contigs with BLASTN (per assembler)
rule blastn:
    input:
        contigs=rules.wrangle_pilon.output.csv,
    output:
        blast=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/BLASTN/{{sample}}.contigs.blastn.csv",
    params:
        temp_file=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/BLASTN/temp.blastn.fasta",
        db=BLASTN_DB,
    threads: THREADS
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/blastn_{{assembler}}.log",
    conda:
        "../envs/blastn.yaml"
    script:
        "../scripts/run_blastn.py"


# Rule: Assess host contamination of contigs with CheckV (per assembler)
rule checkv:
    input:
        contigs=rules.pilon.output.improved_contigs,
    output:
        checkv=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/CHECKV/{{sample}}.contamination.tsv",
    params:
        db=CHECKV_DB,
        folder=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/CHECKV",
    threads: THREADS
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/checkv_{{assembler}}.log",
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
        #
        # A second macOS issue: external volumes formatted as HFS+ or APFS
        # write AppleDouble '._*.hmm' metadata files alongside every real
        # HMM. CheckV mis-reads these as HMMs and fails with "N hmmsearch
        # tasks failed". They are silently removed here on Darwin before
        # CheckV runs; on Linux the block is a no-op.
        """
        if [ "$(uname)" = "Darwin" ]; then
            find {params.db} -name "._*" -delete 2>/dev/null || true
        fi

        checkv contamination \
            -d {params.db} \
            {input.contigs} \
            {params.folder} \
            -t {threads} \
            2> {log}

        mv {params.folder}/contamination.tsv {output.checkv}
        # Drop CheckV intermediates; keep the contamination TSV only.
        # Log a warning rather than failing if any file cannot be removed.
        ls -d -1 {params.folder}/* 2>/dev/null | grep -v .tsv \
            | xargs -r rm -rf \
            || echo "[assembly] WARNING: CheckV cleanup incomplete in {params.folder}" >&2
        """


# Rule: Merge CheckV contamination calls into the BLASTN annotation table
rule merge_checkv_blastn:
    input:
        checkv=rules.checkv.output.checkv,
        blastn=rules.blastn.output.blast,
    output:
        merged_csv=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/CHECKV/{{sample}}.merged.csv",
    params:
        virus_parquet=config["VIRUS_PARQUET"],
        nodes_dmp=config.get("TAXDUMP_NODES", "") or "",
        # `names.dmp` lives next to `nodes.dmp` for parquet refreshes
        # produced by `refresh/refresh_virus_parquet.smk`. Derive it
        # from the nodes path so we only need one config entry.
        names_dmp=(
            (Path(config.get("TAXDUMP_NODES", "")).parent / "names.dmp").as_posix()
            if config.get("TAXDUMP_NODES")
            else ""
        ),
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/panel.yaml"
    run:
        import pandas as pd

        from scripts.functions import canonicalise_blast_match_name

        blastn_df = pd.read_csv(input.blastn)
        checkv_df = (
            pd.read_csv(input.checkv, sep="\t")
            .rename(columns={"contig_id": "name"})
            [["name", "total_genes", "viral_genes", "host_genes", "provirus"]]
        )
        # Inner join on `name` matches the original virusHanter behaviour:
        # contigs without a CheckV entry are dropped from the merged table.
        merged = pd.merge(blastn_df, checkv_df, on="name", how="inner")
        # `assembler` is added by wrangle_pilon and flows through BLASTN;
        # belt-and-braces in case the CSV was rewritten without it.
        if "assembler" not in merged.columns:
            merged = merged.assign(assembler=wildcards.assembler)

        # Canonicalise the BLAST subject title via the parent walk
        # through NCBI's taxdump so the Assembly classification chart
        # no longer renders two bars for what is biologically the
        # same species (the EBV-1 / EBV-2 case, the HSV-1 strain
        # entries, the HHV-6A / HHV-6B split, ...). The function
        # degrades to a no-op when TAXDUMP_NODES is empty or the
        # dmp files are unreadable; `match_name_raw` is always added
        # so the audit trail survives.
        # canonicalise_blast_match_name only consumes the parquet's
        # `name` and `tax_id` columns (via parquet_accession_to_taxid).
        # Project those two so pyarrow skips the multi-hundred-MB
        # `sequence` column on disk; this rule runs once per
        # (sample, assembler), so the saving multiplies. Fall back to a
        # full read if an older parquet build lacks a projected column.
        try:
            try:
                parquet_df = pd.read_parquet(
                    params.virus_parquet, columns=["name", "tax_id"]
                )
            except (ValueError, KeyError):
                parquet_df = pd.read_parquet(params.virus_parquet)
        except Exception as e:
            print(f"[merge_checkv_blastn] could not read parquet: {e}; skipping canonicalisation")
            parquet_df = pd.DataFrame()
        merged = canonicalise_blast_match_name(
            merged,
            parquet_df,
            params.nodes_dmp or None,
            params.names_dmp or None,
        )

        merged.to_csv(output.merged_csv, index=False)


# Rule: Optional second viral-contig classifier (geNomad), per assembler.
#
# Off by default. Turn on with `GENOMAD: "TRUE"` and a populated
# `GENOMAD_DB` in config. One geNomad run per (sample, assembler);
# the per-assembler summary TSV is the headline output:
# `<sample>/<assembler>/GENOMAD/<sample>_summary/<sample>_virus_summary.tsv`.
rule genomad:
    input:
        contigs=rules.pilon.output.improved_contigs,
    output:
        # geNomad names every output directory and file after the
        # input FASTA's stem, not after the sample name. Pilon's
        # improved-contigs output is `<sample>_improved_contigs.fasta`,
        # so geNomad writes to `<sample>_improved_contigs_summary/`.
        # Declaring the rule output to match the actual file avoids
        # the silent "Missing output files" failure that would
        # otherwise fire after geNomad finished successfully.
        summary=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/GENOMAD/{{sample}}_improved_contigs_summary/{{sample}}_improved_contigs_virus_summary.tsv",
    params:
        db=GENOMAD_DB,
        out_dir=f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/GENOMAD",
        # geNomad's `annotate` step calls `mmseqs prefilter` which
        # allocates large amounts of memory in proportion to the
        # query proteome size. metaSPAdes typically produces ~3x
        # more contigs than MEGAHIT, which on a memory-tight host
        # (e.g. an 18 GB laptop) is enough to trigger an OOM SIGKILL.
        # `--splits N` is geNomad's documented memory mitigation:
        # the mmseqs search is partitioned into N chunks and each
        # chunk's peak memory shrinks roughly linearly. Default to
        # 4 splits, which keeps the peak under ~6 GB on the DRRKK
        # samples and still completes in reasonable time. Set to 0
        # in config[GENOMAD_SPLITS] to restore mmseqs' auto-split
        # behaviour on Linux machines with abundant RAM.
        splits=_GENOMAD_SPLITS,
    threads: THREADS
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/genomad_{{assembler}}.log"
    conda:
        "../envs/genomad.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        genomad end-to-end \
            --threads {threads} \
            --splits {params.splits} \
            --cleanup \
            {input.contigs} \
            {params.out_dir} \
            {params.db} \
            > {log} 2>&1
        """
