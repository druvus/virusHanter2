# classification.smk

# Import custom functions
from scripts.functions import (
    canonicalise_taxon_names,
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
        r1=host_removed_r1,
        r2=host_removed_r2,
        fmi=_kaiju_fmi,
        nodes=_kaiju_nodes,
    output:
        kaiju_out=f"{RESULT_FOLDER}/{{sample}}/KAIJU/{{sample}}.kaiju.out",
    threads: THREADS
    resources:
        mem_mb=32000,
        runtime=240,
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
    params:
        # Optional - when TAXDUMP_NODES is set (typically pointing at
        # the refresh-workflow's published nodes.dmp + sibling
        # names.dmp), the kaiju table's ``taxon_name`` column is
        # post-rewritten to the ICTV-binomial species name via a
        # parent-rank walk-up. Degrades to a no-op when unset.
        taxdump_nodes=TAXDUMP_NODES,
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/kaiju.yaml"
    script:
        "../scripts/run_kaiju_to_table.py"

# Rule: Kraken2 classification
rule kraken:
    input:
        r1=host_removed_r1,
        r2=host_removed_r2,
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
    params:
        taxdump_nodes=TAXDUMP_NODES,
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/panel.yaml"
    run:
        df = wrangle_kraken(input.kraken_report)
        # Same species-rank walkup as kaiju_to_table and the BLAST
        # canonicaliser, so every classifier's output uses ICTV-binomial
        # species names. Degrades to a no-op when TAXDUMP_NODES is unset.
        from pathlib import Path as _P
        if params.taxdump_nodes and _P(str(params.taxdump_nodes)).is_file():
            _names = _P(str(params.taxdump_nodes)).parent / "names.dmp"
            if _names.is_file():
                df = canonicalise_taxon_names(
                    df,
                    taxid_col="taxonomy_id",
                    name_col="name",
                    nodes_dmp=str(params.taxdump_nodes),
                    names_dmp=str(_names),
                )
        df.to_csv(output.kraken_csv, index=False)