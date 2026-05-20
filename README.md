# virusHanter2

A Snakemake pipeline for viral metagenomics analysis. Paired-end
Illumina reads (typically Twist Comprehensive Virus Research Panel
enrichment) are quality-trimmed, host-cleaned, classified
(Kraken2 + Kaiju), assembled (MEGAHIT + Pilon), annotated
(BLASTN, CheckV, optionally geNomad), and rendered into an
interactive HTML per sample plus tabular per-batch and
per-(sample, virus) summaries.

`virusHanter2` is the modular refactor of the original `virusHanter`
monolith; HTML rendering is delegated to the
[`reportHanter`](../reportHanter) package.

## Quick start (Linux)

```bash
# 1. driver environment — only snakemake itself needs to be on PATH
#    before the first run; every per-rule tool gets a fresh conda env
#    materialised under .snakemake/conda/ on first use.
conda create -n virushanter -c conda-forge -c bioconda \
    'snakemake-minimal=9.14.*' mamba
conda activate virushanter

# 2. configure
cp config/config.yaml config/config.local.yaml   # then edit DB paths

# 3. dry-run (DAG check only; no tools invoked)
snakemake -n --sdm conda --configfile config/config.local.yaml

# 4. run
snakemake --sdm conda --cores 8 --configfile config/config.local.yaml
```

`snakemake-minimal` is pinned to match what the `reporthanter` rule
env carries, so Snakemake's `script:` directive does not hit
pickle-version mismatches across the driver / per-rule envs.

## Full-feature run on a Linux server

A production run with every optional stage enabled — duplicate-aware
host removal, QUAST assembly QC, geNomad as a second viral-contig
classifier, MultiQC across the batch — looks like this:

```bash
# Driver env
conda create -n virushanter -c conda-forge -c bioconda \
    'snakemake-minimal=9.14.*' mamba
conda activate virushanter

# All reference databases must exist on disk before the run.
# See docs/REFERENCE_DBS.md for sources and refresh cadences.

cat > config/config.prod.yaml <<'YAML'
SAMPLES:        "/data/runs/<run_id>"
RESULTS_FOLDER: "/data/results"
THREADS: 16

HUMAN_INDEX:   "/refs/bwa/human_gencode"
KAIJU_DB:      "/refs/kaiju/refseq"
KRAKEN_DB:     "/refs/kraken2/pluspf"
BLASTN_DB:     "/refs/blast/viral_rna_mito"
CHECKV_DB:     "/refs/checkv/checkv-db-v1.5"
VIRUS_PARQUET: "/refs/individual_virus_fasta/all_viruses.parquet"
GENOMAD_DB:    "/refs/genomad/genomad_db"

CONTIG_LENGTH:  500
NUMBER_OF_PLOTS: 10
COVERAGE_WINDOW: 100
PILON_MEM:      "16G"
MEGAHIT_MEM_FRACTION: 0.8

# Every optional stage turned on.
MULTIQC:     "TRUE"
DEDUPLICATE: "TRUE"
QUAST:       "TRUE"
GENOMAD:     "TRUE"
YAML

snakemake -n --sdm conda --configfile config/config.prod.yaml      # dry-run
snakemake    --sdm conda --cores 16 --configfile config/config.prod.yaml
```

First-run conda env materialisation is the slow step; subsequent runs
reuse the cached envs under `.snakemake/conda/<hash>/`. Add
`--rerun-triggers mtime` if you want only file-timestamp-based
re-runs.

## Documentation

| Topic | File |
|---|---|
| Pipeline stages and output tree | [docs/PIPELINE.md](docs/PIPELINE.md) |
| Config schema, env list, optional flags | [docs/CONFIGURATION.md](docs/CONFIGURATION.md) |
| Reference databases (sources, refresh, rebuild recipes) | [docs/REFERENCE_DBS.md](docs/REFERENCE_DBS.md) |
| Per-(sample, virus) CSV schema + multi-run merge | [docs/PER_VIRUS_OUTPUT.md](docs/PER_VIRUS_OUTPUT.md) |
| Parity invariants with the original `virusHanter` | [docs/PARITY_NOTES.md](docs/PARITY_NOTES.md) |
| Local smoke testing | [test/README.md](test/README.md) |
| Project conventions for AI assistants | [CLAUDE.md](CLAUDE.md) |

## Combining multiple runs

```bash
python scripts/merge_runs.py \
    --result-folder /path/to/RESULTS/<batch1> \
    --result-folder /path/to/RESULTS/<batch2> \
    --out-dir /path/to/master/
# writes master_per_sample.csv + master_per_virus.csv
```

## Support and licence

Open an issue on the repository for problems or questions.
Licensed under the MIT Licence — see [LICENSE](LICENSE).
