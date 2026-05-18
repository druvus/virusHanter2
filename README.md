# virusHanter2

A Snakemake pipeline for viral metagenomics analysis. Paired-end
Illumina reads (typically Twist Comprehensive Virus Research Panel
enrichment) are quality-trimmed, host-cleaned, classified
(Kraken2 + Kaiju), assembled (MEGAHIT + Pilon), annotated
(BLASTN, CheckV, optionally geNomad), and rendered into interactive
HTML per sample plus tabular per-batch and per-(sample, virus)
summaries.

`virusHanter2` is the modular refactor of the original `virusHanter`
monolith; HTML rendering is delegated to the
[`reportHanter`](../reportHanter) package.

## Quick start

```bash
# 1. driver environment
conda create -n virushanter -c bioconda -c conda-forge snakemake
conda activate virushanter

# 2. configure
cp config/config.yaml config/config.local.yaml   # then edit paths

# 3. dry-run (DAG check only)
snakemake -n --sdm conda --configfile config/config.local.yaml

# 4. run
snakemake --sdm conda --cores 8 --configfile config/config.local.yaml
```

Per-rule Conda environments are materialised on first use under
`.snakemake/conda/`.

## Documentation

| Topic | File |
|---|---|
| Config schema, env list, optional flags | [docs/CONFIGURATION.md](docs/CONFIGURATION.md) |
| Pipeline stages and output tree | [docs/PIPELINE.md](docs/PIPELINE.md) |
| Per-(sample, virus) CSV schema + multi-run merge | [docs/PER_VIRUS_OUTPUT.md](docs/PER_VIRUS_OUTPUT.md) |
| Reference databases (sources, refresh, Apple-Silicon caveats) | [docs/REFERENCE_DBS.md](docs/REFERENCE_DBS.md) |
| Parity invariants with the original `virusHanter` | [docs/PARITY_NOTES.md](docs/PARITY_NOTES.md) |
| Local smoke testing | [test/README.md](test/README.md) |

## Combining multiple runs

```bash
python scripts/merge_runs.py \
    --result-folder /path/to/RESULTS/<batch1> \
    --result-folder /path/to/RESULTS/<batch2> \
    --out-dir /path/to/master/
# writes master_per_sample.csv + master_per_virus.csv
```

## Support and license

Open an issue on the repository for problems or questions.
Licensed under the MIT License — see [LICENSE](LICENSE).
