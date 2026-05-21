# refresh/refresh_virus_parquet.smk
#
# Standalone Snakemake workflow that rebuilds VIRUS_PARQUET from
# NCBI Virus. Kept separate from the main pipeline DAG because the
# refresh is an operator task that runs on a different cadence
# (quarterly to half-yearly).
#
# Usage:
#
#   conda activate virushanter
#   cd virusHanter2
#   snakemake -s refresh/refresh_virus_parquet.smk \
#       --configfile refresh/config.yaml --cores 4 \
#       --sdm conda
#
# Re-running with the same config is idempotent thanks to
# Snakemake's mtime tracking. Force a full rebuild with
# `--forcerun build_parquet`.

import os
from pathlib import Path

configfile: "refresh/config.yaml"

OUTPUT_PARQUET = Path(config["OUTPUT_PARQUET"])
DOWNLOAD_DIR = Path(config["DOWNLOAD_DIR"])

# NCBI source URLs. The viral-nucleotide FASTA at
# `viral.1.1.genomic.fna.gz` (and any sibling viral.<N>.1.genomic
# files) is the NCBI Virus all-genomes nucleotide dump. The
# accession2taxid file is needed regardless of source.
NCBI_VIRUS_FASTA_URL = config.get(
    "NCBI_VIRUS_FASTA_URL",
    "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz",
)
ACCESSION2TAXID_URL = config.get(
    "ACCESSION2TAXID_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
)
ACCESSION2TAXID_MD5_URL = ACCESSION2TAXID_URL + ".md5"
TAXDUMP_URL = config.get(
    "TAXDUMP_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
)
TAXDUMP_MD5_URL = TAXDUMP_URL + ".md5"

FASTA_GZ = DOWNLOAD_DIR / "ncbi_virus_nucleotide.fna.gz"
FASTA = DOWNLOAD_DIR / "ncbi_virus_nucleotide.fna"
TAXID_GZ = DOWNLOAD_DIR / "nucl_gb.accession2taxid.gz"
TAXID_MD5 = DOWNLOAD_DIR / "nucl_gb.accession2taxid.gz.md5"
TAXDUMP_TAR = DOWNLOAD_DIR / "taxdump.tar.gz"
TAXDUMP_MD5 = DOWNLOAD_DIR / "taxdump.tar.gz.md5"
TAXDUMP_NODES = DOWNLOAD_DIR / "taxdump" / "nodes.dmp"

# The runtime pipeline points its TAXDUMP_NODES config key at this
# stable location so the rule does not have to re-extract the tar
# every refresh.
TAXDUMP_NODES_PUBLISHED = OUTPUT_PARQUET.parent / "nodes.dmp"


rule all:
    input:
        OUTPUT_PARQUET,
        OUTPUT_PARQUET.with_name(OUTPUT_PARQUET.stem + "_build_stats.json"),
        str(TAXDUMP_NODES_PUBLISHED),


rule download_ncbi_virus_fasta:
    output:
        gz=str(FASTA_GZ),
    params:
        url=NCBI_VIRUS_FASTA_URL,
    log:
        str(DOWNLOAD_DIR / "logs" / "download_ncbi_virus_fasta.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gz}) $(dirname {log})
        wget --quiet --show-progress -O {output.gz} {params.url} 2> {log}
        """


rule decompress_fasta:
    input:
        gz=rules.download_ncbi_virus_fasta.output.gz,
    output:
        fasta=str(FASTA),
    log:
        str(DOWNLOAD_DIR / "logs" / "decompress_fasta.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        gunzip -kc {input.gz} > {output.fasta} 2> {log}
        """


rule download_accession2taxid:
    output:
        gz=str(TAXID_GZ),
        md5=str(TAXID_MD5),
    params:
        url=ACCESSION2TAXID_URL,
        md5_url=ACCESSION2TAXID_MD5_URL,
    log:
        str(DOWNLOAD_DIR / "logs" / "download_accession2taxid.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gz})
        wget --quiet --show-progress -O {output.gz} {params.url} 2> {log}
        wget --quiet --show-progress -O {output.md5} {params.md5_url} 2>> {log}
        # Verify checksum. The .md5 file lists the file as it sits on
        # the NCBI FTP, so md5sum -c only works if the downloaded file
        # has the same basename. cd into the download dir to make that
        # work without rewriting the .md5.
        ( cd $(dirname {output.gz}) && md5sum -c $(basename {output.md5}) ) >> {log} 2>&1
        """


rule download_taxdump:
    output:
        tar=str(TAXDUMP_TAR),
        md5=str(TAXDUMP_MD5),
    params:
        url=TAXDUMP_URL,
        md5_url=TAXDUMP_MD5_URL,
    log:
        str(DOWNLOAD_DIR / "logs" / "download_taxdump.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p $(dirname {output.tar})
        wget --quiet --show-progress -O {output.tar} {params.url} 2> {log}
        wget --quiet --show-progress -O {output.md5} {params.md5_url} 2>> {log}
        ( cd $(dirname {output.tar}) && md5sum -c $(basename {output.md5}) ) >> {log} 2>&1
        """


rule decompress_taxdump:
    input:
        tar=rules.download_taxdump.output.tar,
    output:
        nodes=str(TAXDUMP_NODES),
    params:
        out_dir=str(TAXDUMP_NODES.parent),
    log:
        str(DOWNLOAD_DIR / "logs" / "decompress_taxdump.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        # Extract only nodes.dmp from the tarball; names.dmp and the
        # rest are not needed for rank lookup or genus walk-up.
        tar -xzf {input.tar} -C {params.out_dir} nodes.dmp > {log} 2>&1
        """


rule build_parquet:
    input:
        fasta=rules.decompress_fasta.output.fasta,
        taxid=rules.download_accession2taxid.output.gz,
        taxdump_nodes=rules.decompress_taxdump.output.nodes,
    output:
        parquet=str(OUTPUT_PARQUET),
        stats=str(
            OUTPUT_PARQUET.with_name(OUTPUT_PARQUET.stem + "_build_stats.json")
        ),
    log:
        str(DOWNLOAD_DIR / "logs" / "build_parquet.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        python scripts/build_virus_parquet.py \
            --source ncbi-virus \
            --one-rep-per-taxid \
            --fasta {input.fasta} \
            --taxid {input.taxid} \
            --taxdump-nodes {input.taxdump_nodes} \
            --out {output.parquet} \
            > {log} 2>&1

        # Surface the headline numbers in the Snakemake log.
        echo "--- build_stats.json ---" >> {log}
        cat {output.stats} >> {log}
        """


rule publish_taxdump_nodes:
    """Copy the extracted ``nodes.dmp`` next to the parquet so the
    main pipeline's ``TAXDUMP_NODES`` config key can point at a
    stable location that survives the temporary download workdir."""
    input:
        nodes=rules.decompress_taxdump.output.nodes,
    output:
        nodes=str(TAXDUMP_NODES_PUBLISHED),
    shell:
        """
        mkdir -p $(dirname {output.nodes})
        cp {input.nodes} {output.nodes}
        """
