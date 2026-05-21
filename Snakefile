# Snakefile

# Import necessary modules
import os
from pathlib import Path

# Load configuration
configfile: "config/config.yaml"

# Validate the config against its JSON schema before any rule depends on it.
# Catches the common mistake of running the workflow against the template
# config.yaml without filling in the database paths.
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate
validate(config, "config/config.schema.yaml")


def _check_db_paths_exist(cfg: dict) -> None:
    """Fail fast with a single clean message if any required reference
    database is not on disk where the config claims.

    Catches the common operational mishap of an unmounted external
    drive: without this, the first rule that touches the missing path
    crashes mid-run with a tool-specific error (samtools, blastn,
    etc.) instead of a clear "your database is missing" report. The
    check is best-effort — it tolerates BWA / BLAST prefix-only paths
    by checking the parent directory.
    """
    missing: list[str] = []

    def check_dir(key: str) -> None:
        path = cfg.get(key)
        if not path:
            return
        if not Path(path).is_dir():
            missing.append(f"  {key}: {path}  (expected directory)")

    def check_file(key: str) -> None:
        path = cfg.get(key)
        if not path:
            return
        if not Path(path).is_file():
            missing.append(f"  {key}: {path}  (expected file)")

    def check_prefix(key: str) -> None:
        # BWA / BLAST databases are referenced by a prefix that has no
        # file on disk; the indexed siblings live in the same dir.
        # Treat the parent dir's existence as proxy.
        path = cfg.get(key)
        if not path:
            return
        if not Path(path).parent.is_dir():
            missing.append(
                f"  {key}: {path}  (parent directory missing)"
            )

    check_prefix("HUMAN_INDEX")
    check_dir("KAIJU_DB")
    check_dir("KRAKEN_DB")
    check_prefix("BLASTN_DB")
    check_dir("CHECKV_DB")
    check_file("VIRUS_PARQUET")
    if cfg.get("GENOMAD", "FALSE") == "TRUE":
        check_dir("GENOMAD_DB")
    if cfg.get("SECONDARY_HOST_INDEX"):
        check_prefix("SECONDARY_HOST_INDEX")

    if missing:
        raise WorkflowError(
            "Reference databases are missing or unreachable. "
            "If they live on an external volume, check that it is "
            "mounted before retrying.\n" + "\n".join(missing)
        )


_check_db_paths_exist(config)


# Import pipeline-side helpers. Report rendering and parsing live in the
# reportHanter package; this Snakefile only owns the data-processing side.
from scripts.functions import (
    read_file_as_blob,
    common_suffix,
    paired_reads,
    kaiju_db_files,
    fastx_file_to_df,
    wrangle_kraken,
    run_blastn,
)

# Set sample information
SAMPLES = paired_reads(config["SAMPLES"])
SUFFIX = common_suffix(config["SAMPLES"])

SAMPLES_FOLDER = config["SAMPLES"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(SAMPLES_FOLDER).name)

# Determine if secondary host is specified
SECONDARY_HOST = config.get("SECONDARY_HOST_INDEX", "")
SECONDARY_HOST_NAME = config.get("SECONDARY_HOST_NAME", "")
SECONDARY_HOST_OR_NOT = bool(SECONDARY_HOST)

# Set clean list based on configuration
clean_list = [f"{RESULT_FOLDER}/analysis_done.txt"] if config.get("CLEAN", "FALSE") == "TRUE" else []

# Optional run-level QC. Default on; set MULTIQC: "FALSE" to skip.
RUN_MULTIQC = config.get("MULTIQC", "TRUE") == "TRUE"

# Optional geNomad second viral-contig classifier. Default off so the
# parity invariant holds; flip GENOMAD: "TRUE" and populate GENOMAD_DB
# to opt in.
RUN_GENOMAD_WF = config.get("GENOMAD", "FALSE") == "TRUE"

# Optional QUAST assembly assessment. Default off so the parity
# invariant holds; flip QUAST: "TRUE" to opt in. Output is also fed
# into MultiQC when both are enabled.
RUN_QUAST_WF = config.get("QUAST", "FALSE") == "TRUE"

# Include rule files
include: "rules/pre_processing.smk"
include: "rules/classification.smk"
include: "rules/assembly.smk"
include: "rules/post_processing.smk"

# Trivial rules that should run on the submission host rather than be queued.
# Note: `wrangle_kraken`, `wrangle_pilon`, and `merge_checkv_blastn` need
# pandas/pyfastx from their `conda:` envs, and Snakemake silently ignores
# `conda:` on a localrule. They run as normal jobs.
localrules:
    all,
    clean_everything,

# Define the final targets of the workflow
rule all:
    input:
        expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
        f"{RESULT_FOLDER}/run_information_{Path(SAMPLES_FOLDER).name}.csv",
        # Per-(sample, virus) detail CSV for the collaborator; concatenated
        # across samples by the aggregate_per_virus rule.
        f"{RESULT_FOLDER}/per_virus_{Path(SAMPLES_FOLDER).name}.csv",
        # Per-sample additive QC outputs (do not feed any other rule).
        expand(f"{RESULT_FOLDER}/{{sample}}/logs/human_markdup_stats.txt", sample=SAMPLES),
        expand(f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}.mosdepth.summary.txt", sample=SAMPLES),
        # Run-level QC, gated by MULTIQC config flag (default TRUE).
        [f"{RESULT_FOLDER}/multiqc_report.html"] if RUN_MULTIQC else [],
        # Optional geNomad classifier, gated by GENOMAD config flag.
        (
            expand(
                f"{RESULT_FOLDER}/{{sample}}/GENOMAD/{{sample}}_summary/{{sample}}_virus_summary.tsv",
                sample=SAMPLES,
            )
            if RUN_GENOMAD_WF
            else []
        ),
        # Optional QUAST assembly assessment, gated by QUAST config flag.
        (
            expand(f"{RESULT_FOLDER}/{{sample}}/QUAST/report.tsv", sample=SAMPLES)
            if RUN_QUAST_WF
            else []
        ),
        clean_list,