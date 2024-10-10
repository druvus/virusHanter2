# Snakefile

# Import necessary modules
import os
from pathlib import Path

# Load configuration
configfile: "config/config.yaml"

# Import custom functions
from scripts.functions import (
    read_file_as_blob,
    common_suffix,
    paired_reads,
    kaiju_db_files,
    fastx_file_to_df,
    wrangle_kraken,
    run_blastn,
    parse_bwa_flagstat,
    parse_fastp,
    plot_flagstat,
    plot_kaiju,
    kraken_df,
    plot_kraken,
    plot_blastn,
    alignment_stats,
    panel_report,
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

# Include rule files
include: "rules/pre_processing.smk"
include: "rules/classification.smk"
include: "rules/assembly.smk"
include: "rules/post_processing.smk"

# Define the final targets of the workflow
rule all:
    input:
        expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
        f"{RESULT_FOLDER}/run_information_{Path(SAMPLES_FOLDER).name}.csv",
        clean_list,