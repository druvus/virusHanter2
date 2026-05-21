# refresh/refresh_virus_parquet.smk
#
# Standalone Snakemake workflow that rebuilds VIRUS_PARQUET from
# NCBI viral RefSeq and additionally builds a locally-tuned Kaiju
# FMI index plus an overlap-with-Kraken2 sidecar. Kept separate
# from the main pipeline DAG because the refresh is an operator
# task that runs on a different cadence (quarterly to
# half-yearly).
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
# `--forcerun build_parquet` (or any other rule).

import os
from pathlib import Path

configfile: "refresh/config.yaml"

OUTPUT_PARQUET = Path(config["OUTPUT_PARQUET"])
DOWNLOAD_DIR = Path(config["DOWNLOAD_DIR"])

# Recursive wget against the NCBI RefSeq viral FTP directory pulls
# every viral.N.1.genomic.fna.gz part. Previous builds only pulled
# the .1.1 part by accident, which dropped roughly half the
# release. The same pattern applies to the protein FASTAs used for
# the Kaiju build (`viral.N.protein.faa.gz`).
REFSEQ_VIRAL_URL_DIR = config.get(
    "REFSEQ_VIRAL_URL_DIR",
    "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/",
)
ACCESSION2TAXID_URL = config.get(
    "ACCESSION2TAXID_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
)
ACCESSION2TAXID_MD5_URL = ACCESSION2TAXID_URL + ".md5"
PROT_ACCESSION2TAXID_URL = config.get(
    "PROT_ACCESSION2TAXID_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz",
)
PROT_ACCESSION2TAXID_MD5_URL = PROT_ACCESSION2TAXID_URL + ".md5"
TAXDUMP_URL = config.get(
    "TAXDUMP_URL",
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
)
TAXDUMP_MD5_URL = TAXDUMP_URL + ".md5"

# Where the production Kraken2 viral DB lives. The overlap diagnostic
# inspects this DB to compare its taxid universe against the parquet's.
KRAKEN_DB_FOR_COMPARE = config.get("KRAKEN_DB_FOR_COMPARE", "")

FASTA = DOWNLOAD_DIR / "ncbi_virus_nucleotide.fna"
FASTA_PARTS_DIR = DOWNLOAD_DIR / "refseq_viral_parts"
FASTA_MARKER = DOWNLOAD_DIR / ".refseq_viral_downloaded"
PROTEIN_FAA = DOWNLOAD_DIR / "refseq_viral_proteins.faa"
PROTEIN_PARTS_DIR = DOWNLOAD_DIR / "refseq_viral_proteins"
PROTEIN_MARKER = DOWNLOAD_DIR / ".refseq_viral_proteins_downloaded"
TAXID_GZ = DOWNLOAD_DIR / "nucl_gb.accession2taxid.gz"
TAXID_MD5 = DOWNLOAD_DIR / "nucl_gb.accession2taxid.gz.md5"
PROT_TAXID_GZ = DOWNLOAD_DIR / "prot.accession2taxid.gz"
PROT_TAXID_MD5 = DOWNLOAD_DIR / "prot.accession2taxid.gz.md5"
TAXDUMP_TAR = DOWNLOAD_DIR / "taxdump.tar.gz"
TAXDUMP_MD5 = DOWNLOAD_DIR / "taxdump.tar.gz.md5"
TAXDUMP_NODES = DOWNLOAD_DIR / "taxdump" / "nodes.dmp"
TAXDUMP_NAMES = DOWNLOAD_DIR / "taxdump" / "names.dmp"

# Kaiju build artefacts.
KAIJU_BUILD_PREFIX = DOWNLOAD_DIR / "kaiju_refseq_viral"
KAIJU_BUILD_FMI = DOWNLOAD_DIR / "kaiju_refseq_viral.fmi"
KAIJU_BUILD_FAA = DOWNLOAD_DIR / "kaiju_refseq_viral.faa"

# Published paths next to the parquet so the main pipeline's
# TAXDUMP_NODES and KAIJU_DB config keys can point at stable
# locations that survive the temporary download workdir.
PARQUET_PARENT = OUTPUT_PARQUET.parent
TAXDUMP_NODES_PUBLISHED = PARQUET_PARENT / "nodes.dmp"
KAIJU_PUBLISH_DIR = PARQUET_PARENT / "kaiju_refseq_viral"
KAIJU_PUBLISHED_FMI = KAIJU_PUBLISH_DIR / "kaiju_refseq_viral.fmi"
KAIJU_PUBLISHED_NODES = KAIJU_PUBLISH_DIR / "nodes.dmp"
KAIJU_PUBLISHED_NAMES = KAIJU_PUBLISH_DIR / "names.dmp"

# Overlap-with-Kraken2 sidecar lives next to the parquet.
OVERLAP_TSV = OUTPUT_PARQUET.with_name(OUTPUT_PARQUET.stem + "_vs_kraken2.tsv")


rule all:
    input:
        OUTPUT_PARQUET,
        OUTPUT_PARQUET.with_name(OUTPUT_PARQUET.stem + "_build_stats.json"),
        str(TAXDUMP_NODES_PUBLISHED),
        str(KAIJU_PUBLISHED_FMI),
        str(KAIJU_PUBLISHED_NODES),
        str(KAIJU_PUBLISHED_NAMES),
        str(OVERLAP_TSV),


rule download_refseq_viral_nucleotide:
    """Pull every viral.N.1.genomic.fna.gz part from the NCBI
    RefSeq viral FTP directory. NCBI publishes the release as
    multiple files; the part count varies per release and cannot
    be predicted at workflow-build time, so the rule fetches the
    directory listing, greps for matching filenames, and downloads
    each. Recursive wget with `--accept '<glob>'` failed silently
    on this NCBI listing (only index.html + robots.txt got pulled),
    hence the enumerate-then-fetch approach."""
    output:
        marker=str(FASTA_MARKER),
    params:
        url_dir=REFSEQ_VIRAL_URL_DIR,
        out_dir=str(FASTA_PARTS_DIR),
    log:
        str(DOWNLOAD_DIR / "logs" / "download_refseq_viral_nucleotide.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p {params.out_dir} $(dirname {log})
        # Fetch the directory listing, extract matching filenames,
        # download each part. `grep -oE` finds the bare filename
        # tokens regardless of the surrounding HTML.
        curl --silent --show-error --fail {params.url_dir} \
            | grep -oE 'viral\\.[0-9]+\\.1\\.genomic\\.fna\\.gz' \
            | sort -u \
            | tee {params.out_dir}/.parts_list > {log}
        if [ ! -s {params.out_dir}/.parts_list ]; then
            echo "ERROR: no viral.*.1.genomic.fna.gz parts found at {params.url_dir}" >> {log}
            exit 1
        fi
        while read part; do
            echo "fetching $part" >> {log}
            wget --no-verbose -O {params.out_dir}/$part \
                {params.url_dir}$part >> {log} 2>&1
        done < {params.out_dir}/.parts_list
        touch {output.marker}
        """


rule decompress_fasta:
    input:
        marker=rules.download_refseq_viral_nucleotide.output.marker,
    output:
        fasta=str(FASTA),
    params:
        in_dir=str(FASTA_PARTS_DIR),
    log:
        str(DOWNLOAD_DIR / "logs" / "decompress_fasta.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        # Concatenate every viral.*.1.genomic.fna.gz part into a
        # single FASTA the builder can stream. zcat is fine on the
        # multi-file glob because the parts are independent
        # gzipped streams.
        zcat {params.in_dir}/viral.*.1.genomic.fna.gz > {output.fasta} 2> {log}
        """


rule download_refseq_viral_proteins:
    """Pull every viral.N.protein.faa.gz part from the NCBI RefSeq
    viral FTP directory. Used to build the matching Kaiju index
    from the same RefSeq snapshot as the parquet. Same
    enumerate-then-fetch approach as
    `download_refseq_viral_nucleotide`."""
    output:
        marker=str(PROTEIN_MARKER),
    params:
        url_dir=REFSEQ_VIRAL_URL_DIR,
        out_dir=str(PROTEIN_PARTS_DIR),
    log:
        str(DOWNLOAD_DIR / "logs" / "download_refseq_viral_proteins.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        curl --silent --show-error --fail {params.url_dir} \
            | grep -oE 'viral\\.[0-9]+\\.protein\\.faa\\.gz' \
            | sort -u \
            | tee {params.out_dir}/.parts_list > {log}
        if [ ! -s {params.out_dir}/.parts_list ]; then
            echo "ERROR: no viral.*.protein.faa.gz parts found at {params.url_dir}" >> {log}
            exit 1
        fi
        while read part; do
            echo "fetching $part" >> {log}
            wget --no-verbose -O {params.out_dir}/$part \
                {params.url_dir}$part >> {log} 2>&1
        done < {params.out_dir}/.parts_list
        touch {output.marker}
        """


rule decompress_viral_proteins:
    input:
        marker=rules.download_refseq_viral_proteins.output.marker,
    output:
        faa=str(PROTEIN_FAA),
    params:
        in_dir=str(PROTEIN_PARTS_DIR),
    log:
        str(DOWNLOAD_DIR / "logs" / "decompress_viral_proteins.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        zcat {params.in_dir}/viral.*.protein.faa.gz > {output.faa} 2> {log}
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
        wget --no-verbose -O {output.gz} {params.url} 2> {log}
        wget --no-verbose -O {output.md5} {params.md5_url} 2>> {log}
        # Verify checksum. The .md5 file lists the file as it sits on
        # the NCBI FTP, so md5sum -c only works if the downloaded file
        # has the same basename. cd into the download dir to make that
        # work without rewriting the .md5.
        ( cd $(dirname {output.gz}) && md5sum -c $(basename {output.md5}) ) >> {log} 2>&1
        """


rule download_prot_accession2taxid:
    """Protein-accession-to-taxid mapping. Drives the header
    rewrite step that prepares the RefSeq viral protein FASTA for
    Kaiju's BWT builder."""
    output:
        gz=str(PROT_TAXID_GZ),
        md5=str(PROT_TAXID_MD5),
    params:
        url=PROT_ACCESSION2TAXID_URL,
        md5_url=PROT_ACCESSION2TAXID_MD5_URL,
    log:
        str(DOWNLOAD_DIR / "logs" / "download_prot_accession2taxid.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gz})
        wget --no-verbose -O {output.gz} {params.url} 2> {log}
        wget --no-verbose -O {output.md5} {params.md5_url} 2>> {log}
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
        wget --no-verbose -O {output.tar} {params.url} 2> {log}
        wget --no-verbose -O {output.md5} {params.md5_url} 2>> {log}
        ( cd $(dirname {output.tar}) && md5sum -c $(basename {output.md5}) ) >> {log} 2>&1
        """


rule decompress_taxdump:
    input:
        tar=rules.download_taxdump.output.tar,
    output:
        nodes=str(TAXDUMP_NODES),
        names=str(TAXDUMP_NAMES),
    params:
        out_dir=str(TAXDUMP_NODES.parent),
    log:
        str(DOWNLOAD_DIR / "logs" / "decompress_taxdump.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        # Extract nodes.dmp (for rank / genus walk-up) and names.dmp
        # (Kaiju needs it at classification time alongside nodes.dmp).
        tar -xzf {input.tar} -C {params.out_dir} nodes.dmp names.dmp > {log} 2>&1
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
            --source refseq \
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


rule build_kaiju_refseq_viral:
    """Build a Kaiju FMI index from the concatenated RefSeq viral
    protein FASTAs. The header-rewriter script reformats each
    record's header to Kaiju's expected ``>kaiju|<taxid>|<accession>``
    shape using the protein accession2taxid mapping."""
    input:
        faa=rules.decompress_viral_proteins.output.faa,
        prot_taxid=rules.download_prot_accession2taxid.output.gz,
    output:
        fmi=str(KAIJU_BUILD_FMI),
        reformatted_faa=str(KAIJU_BUILD_FAA),
    params:
        prefix=str(KAIJU_BUILD_PREFIX),
    threads: 4
    log:
        str(DOWNLOAD_DIR / "logs" / "build_kaiju.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        python scripts/reformat_kaiju_headers.py \
            --proteins {input.faa} \
            --prot-taxid {input.prot_taxid} \
            --out {output.reformatted_faa} \
            > {log} 2>&1

        # kaiju-mkbwt builds the Burrows-Wheeler transform; -n
        # thread count, -o output prefix, last argument is the
        # input FASTA. Produces <prefix>.bwt and <prefix>.sa.
        kaiju-mkbwt -n {threads} -o {params.prefix} \
            {output.reformatted_faa} >> {log} 2>&1

        # kaiju-mkfmi consumes the BWT and emits the FM-index.
        kaiju-mkfmi {params.prefix} >> {log} 2>&1
        """


rule publish_kaiju_refseq_viral:
    """Copy the built Kaiju FMI and the taxdump dmps next to the
    parquet so the main pipeline can point its KAIJU_DB config
    key at a stable location."""
    input:
        fmi=rules.build_kaiju_refseq_viral.output.fmi,
        nodes=rules.decompress_taxdump.output.nodes,
        names=rules.decompress_taxdump.output.names,
    output:
        fmi=str(KAIJU_PUBLISHED_FMI),
        nodes=str(KAIJU_PUBLISHED_NODES),
        names=str(KAIJU_PUBLISHED_NAMES),
    shell:
        """
        mkdir -p $(dirname {output.fmi})
        cp {input.fmi} {output.fmi}
        cp {input.nodes} {output.nodes}
        cp {input.names} {output.names}
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


rule compare_with_kraken2:
    """Compute the overlap between the parquet's taxid universe and
    the production Kraken2 viral DB. Emits a sidecar TSV (one row
    per taxid in either set) and extends build_stats.json with
    intersection / parquet-only / kraken2-only counters."""
    input:
        parquet=rules.build_parquet.output.parquet,
        stats=rules.build_parquet.output.stats,
    output:
        overlap_tsv=str(OVERLAP_TSV),
    params:
        kraken_db=KRAKEN_DB_FOR_COMPARE,
    log:
        str(DOWNLOAD_DIR / "logs" / "compare_with_kraken2.log"),
    conda:
        "../envs/refresh.yaml"
    shell:
        """
        python scripts/compare_parquet_kraken2.py \
            --parquet {input.parquet} \
            --kraken-db {params.kraken_db} \
            --build-stats {input.stats} \
            --out {output.overlap_tsv} \
            > {log} 2>&1
        """
