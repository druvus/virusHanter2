# test/build_fixtures.smk
#
# Standalone Snakemake workflow that materialises every mini reference
# database the smoke pipeline consumes. Each rule reuses the same conda env
# YAMLs as the production pipeline, so it picks up the same tool versions.
#
# Invoke from the repo root:
#
#   snakemake --snakefile test/build_fixtures.smk --sdm conda --cores 2 \
#             --directory .
#
# Outputs:
#   test/test_R{1,2}.fastq.gz              (paired reads, gzipped to match
#                                           production .fastq.gz inputs)
#   test/mini_db/host.fasta                (host reference)
#   test/mini_db/virus.fasta               (viral reference, nucleotide)
#   test/mini_db/virus_aa.fasta            (viral reference, protein)
#   test/mini_db/human/human.*             (BWA index files)
#   test/mini_db/kraken/{hash,opts,taxo}.k2d
#   test/mini_db/kaiju/kaiju_db.fmi + nodes.dmp + names.dmp
#   test/mini_db/blast/viral.{nhr,nin,nsq,...}
#   test/mini_db/virus.parquet
#   test/mini_db/checkv/.stub              (placeholder; CheckV not built)

import os

MINI = "test/mini_db"
SCRIPTS = "test/scripts"


rule all:
    input:
        "test/test_R1.fastq.gz",
        "test/test_R2.fastq.gz",
        f"{MINI}/human/human.bwt",
        f"{MINI}/kraken/hash.k2d",
        f"{MINI}/kaiju/kaiju_db.fmi",
        f"{MINI}/blast/viral.nhr",
        f"{MINI}/virus.parquet",
    # CheckV's database is too large to synthesise. The smoke runner
    # detects a real CheckV DB at `test/mini_db/checkv` by looking for
    # `genome_db/`; when absent it degrades to a BLASTN-only smoke.


# Pipeline-side helpers (Python only — no bioinformatics tools required).
rule synthesize_reads:
    output:
        r1="test/test_R1.fastq.gz",
        r2="test/test_R2.fastq.gz",
    shell:
        "python {SCRIPTS}/synthesize_fastq.py "
        "--out-r1 {output.r1} --out-r2 {output.r2}"


rule write_references:
    output:
        host=f"{MINI}/host.fasta",
        virus=f"{MINI}/virus.fasta",
        virus_aa=f"{MINI}/virus_aa.fasta",
    shell:
        "python {SCRIPTS}/write_references.py "
        "--host-fasta {output.host} "
        "--virus-fasta {output.virus} "
        "--virus-protein-fasta {output.virus_aa}"


rule make_virus_parquet:
    output:
        parquet=f"{MINI}/virus.parquet",
    conda:
        "../envs/panel.yaml"
    shell:
        "python {SCRIPTS}/make_virus_parquet.py --out {output.parquet}"


# BWA host index.
rule bwa_host_index:
    input:
        host=rules.write_references.output.host,
    output:
        bwt=f"{MINI}/human/human.bwt",
    params:
        prefix=f"{MINI}/human/human",
    conda:
        "../envs/bwa.yaml"
    shell:
        "mkdir -p $(dirname {params.prefix}) && "
        "bwa index -p {params.prefix} {input.host}"


# Kraken2 mini-DB. We hand-write a 3-line nodes.dmp + names.dmp and tag the
# virus reference with the kraken-style accession comment so kraken2-build
# can pull a taxid out of it.
rule kraken_taxonomy:
    output:
        nodes=f"{MINI}/kraken/taxonomy/nodes.dmp",
        names=f"{MINI}/kraken/taxonomy/names.dmp",
    run:
        os.makedirs(os.path.dirname(output.nodes), exist_ok=True)
        with open(output.nodes, "w") as f:
            f.write(
                "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
                "10239\t|\t1\t|\tsuperkingdom\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
                "100001\t|\t10239\t|\tspecies\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            )
        with open(output.names, "w") as f:
            f.write(
                "1\t|\troot\t|\t\t|\tscientific name\t|\n"
                "10239\t|\tViruses\t|\t\t|\tscientific name\t|\n"
                "100001\t|\tsynthetic virus\t|\t\t|\tscientific name\t|\n"
            )


rule kraken_library:
    input:
        virus=rules.write_references.output.virus,
    output:
        fna=f"{MINI}/kraken/library/virus.fna",
    shell:
        r"""mkdir -p $(dirname {output.fna})
        sed 's/^>synthetic_virus/>synthetic_virus|kraken:taxid|100001/' \
            {input.virus} > {output.fna}"""


rule kraken_build:
    input:
        nodes=rules.kraken_taxonomy.output.nodes,
        names=rules.kraken_taxonomy.output.names,
        fna=rules.kraken_library.output.fna,
    output:
        hashdb=f"{MINI}/kraken/hash.k2d",
    params:
        db=f"{MINI}/kraken",
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken2-build --add-to-library {input.fna} --db {params.db} --no-masking && "
        "kraken2-build --build --db {params.db}"


# Kaiju mini-DB. nodes.dmp / names.dmp are reused from the Kraken taxdump.
rule kaiju_dumps:
    input:
        nodes=rules.kraken_taxonomy.output.nodes,
        names=rules.kraken_taxonomy.output.names,
    output:
        nodes=f"{MINI}/kaiju/nodes.dmp",
        names=f"{MINI}/kaiju/names.dmp",
    shell:
        "mkdir -p $(dirname {output.nodes}) && "
        "cp {input.nodes} {output.nodes} && "
        "cp {input.names} {output.names}"


rule kaiju_protein_fasta:
    input:
        aa=rules.write_references.output.virus_aa,
    output:
        faa=f"{MINI}/kaiju/virus_aa.faa",
    shell:
        r"""sed 's/^>synthetic_virus_aa/>P00001_100001/' {input.aa} > {output.faa}"""


rule kaiju_build:
    input:
        faa=rules.kaiju_protein_fasta.output.faa,
        nodes=rules.kaiju_dumps.output.nodes,
        names=rules.kaiju_dumps.output.names,
    output:
        fmi=f"{MINI}/kaiju/kaiju_db.fmi",
    params:
        outdir=f"{MINI}/kaiju",
    conda:
        "../envs/kaiju.yaml"
    shell:
        "cd {params.outdir} && "
        "kaiju-mkbwt -n 2 -o kaiju_db virus_aa.faa && "
        "kaiju-mkfmi kaiju_db && "
        "rm -f kaiju_db.bwt kaiju_db.sa"


# BLAST nt mini-DB.
rule blast_db:
    input:
        virus=rules.write_references.output.virus,
    output:
        nhr=f"{MINI}/blast/viral.nhr",
    params:
        prefix=f"{MINI}/blast/viral",
    conda:
        "../envs/blastn.yaml"
    shell:
        "mkdir -p $(dirname {params.prefix}) && "
        "makeblastdb -in {input.virus} -dbtype nucl -out {params.prefix} -title viral_mini"


