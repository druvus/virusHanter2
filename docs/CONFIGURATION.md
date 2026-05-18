# Configuration

The pipeline reads two files at startup: `config/config.yaml` (your
inputs and database paths) and `config/config.schema.yaml` (a JSON
schema that fails fast on placeholder paths). All Conda environments
are declared in `envs/*.yaml` and materialised on demand by Snakemake.

## `config/config.yaml`

Required keys:

| Key | Purpose |
|---|---|
| `SAMPLES` | Directory containing paired-end FASTQs. Sample names are derived from the common prefix of each `_R1`/`_R2` pair; both `.fastq` and `.fastq.gz` inputs are picked up. |
| `RESULTS_FOLDER` | Output root. The actual results live at `{RESULTS_FOLDER}/{basename(SAMPLES)}/`. |
| `THREADS` | Default per-rule thread count; combined with `snakemake --cores N` for scheduling. |
| `HUMAN_INDEX` | BWA index prefix for the host (human) genome. |
| `KAIJU_DB` | Directory containing `.fmi`, `names.dmp`, `nodes.dmp`. |
| `KRAKEN_DB` | Kraken2 database directory (with `hash.k2d`, `opts.k2d`, `taxo.k2d`). |
| `BLASTN_DB` | BLAST nucleotide database **prefix** (e.g. the `.nal` alias name). The runner derives `BLASTDB` from the parent directory so taxdb lookups work. |
| `CHECKV_DB` | CheckV database directory. |
| `VIRUS_PARQUET` | Parquet with columns `(name, sequence, tax_id)` used to pick references for the kraken-top-N coverage step. See [REFERENCE_DBS.md](REFERENCE_DBS.md) for the build recipe. |
| `PLOT_THRESHOLD` | Minimum mean coverage for a reference to appear in the SVG coverage plots (default `5`). |
| `NUMBER_OF_PLOTS` | Top-N Kraken viral hits to include (default `10`). Same cap drives the per-virus CSV row count. |

Optional keys:

| Key | Default | Purpose |
|---|---|---|
| `CLEAN` | `"FALSE"` | If `"TRUE"`, remove intermediates after the run and write `analysis_done.txt`. |
| `CONTIG_LENGTH` | `500` | Minimum polished-contig length kept after Pilon. |
| `PILON_MEM` | `"50G"` | JVM heap for Pilon. |
| `MEGAHIT_MEM_FRACTION` | `0.5` | Fraction of system RAM MEGAHIT is allowed to allocate. Drop to `0.2-0.3` on memory-tight laptops. |
| `MULTIQC` | `"TRUE"` | Emit `{batch}/multiqc_report.html` at the end of the run. |
| `GENOMAD` | `"FALSE"` | Run geNomad alongside CheckV. Requires `GENOMAD_DB`. |
| `GENOMAD_DB` | `""` | Path to the populated `genomad_db/` directory. See [REFERENCE_DBS.md](REFERENCE_DBS.md#sources-for-the-inputs). |
| `SECONDARY_HOST_INDEX` | unset | BWA prefix for a second host (e.g. mouse). Adds a second host-removal stage when set. |
| `SECONDARY_HOST_NAME` | unset | Display name shown in the per-sample report. |

Use absolute paths for everything. The schema rejects the
placeholder `/path/to/...` strings shipped in `config.yaml.example`
so a misconfigured run fails immediately.

## Production config

`config/config.production.yaml` is the workstation-ready config
pointing at the LaCie reference databases at
`/Volumes/LaCie/REGIONEN/ref_dbs/`. Use it as a template:

```
snakemake --sdm conda --cores 4 --configfile config/config.production.yaml
```

## Conda environments

Each rule declares its own env in `envs/`:

| File | Tools |
|---|---|
| `envs/fastp.yaml` | fastp |
| `envs/bwa.yaml` | bwa, samtools |
| `envs/samtools.yaml` | samtools (used by `markdup_human`, `remove_host`, etc.) |
| `envs/kraken.yaml` | kraken2 |
| `envs/kaiju.yaml` | kaiju |
| `envs/megahit.yaml` | megahit |
| `envs/pilon.yaml` | pilon, bwa, samtools, openjdk |
| `envs/blastn.yaml` | blast, pandas, pyfastx |
| `envs/checkv.yaml` | checkv |
| `envs/bam2plot.yaml` | bam2plot (pip), samtools, polars-lts-cpu |
| `envs/mosdepth.yaml` | mosdepth |
| `envs/multiqc.yaml` | multiqc |
| `envs/genomad.yaml` | genomad |
| `envs/panel.yaml` | python, pandas, pyfastx, pyarrow (wrangling rules) |
| `envs/reporthanter.yaml` | python 3.12 + pip-installed `reportHanter` from GitHub |

Snakemake materialises an env the first time any rule that declares
it runs. Subsequent runs reuse the cached env under
`.snakemake/conda/<hash>/`.

## Apple Silicon notes

The driver env stays native osx-arm64. A few bioconda tools have
historical or current rough edges on this platform; see
[../test/run_smoke.sh](../test/run_smoke.sh) for the explicit list
of workarounds (MEGAHIT thread cap, `--no-hw-accel`, RAM-bound
mem fraction).
