# classification.smk

# Import custom functions
from scripts.functions import (
    kaiju_db_files,
    wrangle_kraken,
)

# Set variables
THREADS = config["THREADS"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(config["SAMPLES"]).name)
KAIJU_DB = config["KAIJU_DB"]
KRAKEN_DB = config["KRAKEN_DB"]


# Lazy input resolvers for the Kaiju database files. These run at DAG build
# time (after `--lint` / `-n` have already succeeded), so a missing database
# directory will not abort workflow parsing.
def _kaiju_fmi(wildcards):
    return str(kaiju_db_files(config["KAIJU_DB"])[0])


def _kaiju_names(wildcards):
    return str(kaiju_db_files(config["KAIJU_DB"])[1])


def _kaiju_nodes(wildcards):
    return str(kaiju_db_files(config["KAIJU_DB"])[2])


# Rule: Kaiju classification
rule kaiju:
    input:
        r1=lambda wildcards: host_removed_r1(wildcards),
        r2=lambda wildcards: host_removed_r2(wildcards),
        fmi=_kaiju_fmi,
        nodes=_kaiju_nodes,
    output:
        kaiju_out=f"{RESULT_FOLDER}/{{sample}}/KAIJU/{{sample}}.kaiju.out",
    threads: THREADS
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/kaiju.log"
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju \
            -t {input.nodes} \
            -f {input.fmi} \
            -i {input.r1} \
            -j {input.r2} \
            -z {threads} \
            -o {output.kaiju_out} \
            > {log} 2>&1
        """

# Rule: Kaiju output to table
rule kaiju_to_table:
    input:
        kaiju_out=rules.kaiju.output.kaiju_out,
        names=_kaiju_names,
        nodes=_kaiju_nodes,
    output:
        kaiju_table=f"{RESULT_FOLDER}/{{sample}}/KAIJU/{{sample}}.kaiju.table.tsv",
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju2table \
            -t {input.nodes} \
            -n {input.names} \
            -r genus \
            -e \
            -o {output.kaiju_table} \
            {input.kaiju_out}
        """

# Rule: Kraken2 classification
rule kraken:
    input:
        r1=lambda wildcards: host_removed_r1(wildcards),
        r2=lambda wildcards: host_removed_r2(wildcards),
    output:
        kraken_report=f"{RESULT_FOLDER}/{{sample}}/KRAKEN/{{sample}}.kraken.report",
    params:
        db=KRAKEN_DB,
    threads: THREADS
    resources:
        mem_mb=64000,
        runtime=120,
    log:
        f"{RESULT_FOLDER}/{{sample}}/logs/kraken.log"
    conda:
        "../envs/kraken.yaml"
    shell:
        """
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --report {output.kraken_report} \
            --use-names \
            --paired \
            {input.r1} \
            {input.r2} \
            > {log} 2>&1
        """

# Rule: Process Kraken2 output
rule wrangle_kraken:
    input:
        kraken_report=rules.kraken.output.kraken_report,
    output:
        kraken_csv=f"{RESULT_FOLDER}/{{sample}}/KRAKEN/{{sample}}.kraken.csv",
    conda:
        "../envs/panel.yaml"
    run:
        df = wrangle_kraken(input.kraken_report)
        df.to_csv(output.kraken_csv, index=False)