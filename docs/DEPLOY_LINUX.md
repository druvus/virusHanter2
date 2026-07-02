# Deploying on a fresh Linux server

End-to-end setup for processing many sequencing runs on a new Linux
host: install, build the reference databases, run the pipeline per run,
and produce a single merged output across all runs.

Companion docs: [DATABASE_SETUP.md](DATABASE_SETUP.md) (per-database
detail), [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md) (the database build
workflow), [CONFIGURATION.md](CONFIGURATION.md) (every config key),
[PARITY_SIGNOFF.md](PARITY_SIGNOFF.md) (validating against the original
virusHanter).

## Set one base folder

Everything below derives from a single environment variable. Choose it
once; export it in your shell profile to make it persist.

```bash
export VH2_ROOT=/data          # <-- the only path you must choose
export REF="$VH2_ROOT/refdbs"  # derived; used throughout
```

Layout used throughout this guide:

```
$VH2_ROOT/
  code/      cloned repositories
  refdbs/    reference databases  (= $REF)
  runs/      one subfolder per sequencing run (your input FASTQs)
  results/   pipeline output, one subfolder per run
  master/    merged cross-run output
```

Create it:

```bash
mkdir -p "$VH2_ROOT"/{code,refdbs,runs,results,master}
```

The example config files ship with `/data` paths. The copy steps below
rewrite `/data` to `$VH2_ROOT` automatically with `sed`, so keeping
`VH2_ROOT=/data` changes nothing and any other base just works. (If you
need the databases and the raw runs on different volumes, edit the
individual paths in the configs by hand instead of relying on the
single base.)

## Architecture note (shapes the run/merge steps)

The pipeline processes **one `SAMPLES` directory (one run) per
invocation**. Sample discovery is a flat, regex scan of that directory,
not recursive. So with N run folders you invoke the pipeline N times
(one per run) and then merge the per-run CSVs with
`scripts/merge_runs.py`.

Each run folder must contain **only** the paired FASTQs for its samples
(`*.fastq.gz` / `*.fq.gz` and similar), an **even** count, with R1/R2
sorting adjacently (e.g. `sampleA_R1_001.fastq.gz`,
`sampleA_R2_001.fastq.gz`). An odd count raises a clear error. Place (or
symlink) your runs under `$VH2_ROOT/runs/`.

## 0. Prerequisites

- Linux (required: CheckV 1.0.3 misbehaves on macOS).
- conda / mamba (miniforge recommended).
- Disk: ~40 GB scratch for the database build, ~5 GB published viral DBs,
  ~14 GB human BWA index, ~2 GB CheckV, plus results. RAM: the
  viral-scoped Kraken2 / Kaiju built below are light (a few GB); the
  human BWA index needs ~4 GB at run time.

## 1. Install

```bash
# Driver/conda env. snakemake-minimal is pinned to 9.23.* for pickle
# compatibility with the reporthanter rule env. The run:-block rules
# (wrangle_pilon, merge_checkv_blastn, wrangle_kraken) execute in THIS
# env, so it also needs pandas / numpy / pyfastx / pyarrow.
mamba create -n virushanter -c conda-forge -c bioconda \
  'snakemake-minimal=9.23.*' 'pandas>=3' 'numpy>=1.24' 'pyfastx>=2.0' 'pyarrow>=14'
conda activate virushanter

cd "$VH2_ROOT/code"
git clone https://github.com/druvus/virusHanter2.git
cd virusHanter2
```

Every per-tool step uses its own conda env from `envs/*.yaml` via
`--sdm conda` (created on first use). reportHanter is installed by its
rule env at the pinned `@v0.9.0` tag; nothing to install by hand.

## 2. Build the NCBI-virus databases (BLAST + Kaiju + Kraken2 + parquet + taxdump)

One coordinated snapshot for all four classifier databases plus
`nodes.dmp` / `names.dmp`, via the refresh workflow. Copying the example
config rewrites `/data` to your `$VH2_ROOT`:

```bash
sed "s|/data|$VH2_ROOT|g" refresh/config.local.example.yaml > refresh/config.local.yaml
$EDITOR refresh/config.local.yaml      # review OUTPUT_PARQUET, DOWNLOAD_DIR, KRAKEN_DB_FOR_COMPARE

snakemake -s refresh/refresh_virus_parquet.smk \
  --configfile refresh/config.local.yaml --cores 4 --sdm conda \
  --omit-from compare_with_kraken2     # skip the optional overlap sidecar on a fresh server
# (optional, later) re-run WITHOUT --omit-from to write all_viruses_vs_kraken2.tsv
```

Published under `$REF/virus_ref/` (these map to the main config keys):

```
$REF/virus_ref/
  all_viruses.parquet                       -> VIRUS_PARQUET
  nodes.dmp  names.dmp                       -> TAXDUMP_NODES (nodes.dmp)
  kaiju_refseq_viral/  (.fmi+nodes+names)    -> KAIJU_DB
  kraken2_refseq_viral/ (hash/taxo/opts...)  -> KRAKEN_DB
  blast_refseq_viral/viral_rna_mito(.nal)    -> BLASTN_DB (prefix, no extension)
```

These Kraken2 / Kaiju DBs are viral-scoped: coordinated and light, ideal
for viral metagenomics. If you later need bacterial/fungal context, point
`KRAKEN_DB` at a pre-built `pluspf` instead (see DATABASE_SETUP.md); not
required for viral detection.

## 3. The two databases the refresh does not build

Human BWA index (`HUMAN_INDEX`):

```bash
mkdir -p "$REF/human" && cd "$REF/human"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mamba create -n bwa -c bioconda -c conda-forge 'bwa>=0.7.17'
conda run -n bwa bwa index -p human_gencode GRCh38.primary_assembly.genome.fa  # ~1.5 h, ~14 GB
# HUMAN_INDEX = $REF/human/human_gencode   (the prefix, not the .fa)
```

CheckV (`CHECKV_DB`):

```bash
mkdir -p "$REF/checkv" && cd "$REF/checkv"
wget https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar -xzf checkv-db-v1.5.tar.gz
# CHECKV_DB = $REF/checkv/checkv-db-v1.5
```

## 4. Configure the main pipeline

```bash
cd "$VH2_ROOT/code/virusHanter2"
sed "s|/data|$VH2_ROOT|g" config/config.local.example.yaml > config/config.local.yaml
$EDITOR config/config.local.yaml       # review DB paths + RESULTS_FOLDER + THREADS
```

The example wires the five refresh-built keys plus `HUMAN_INDEX` and
`CHECKV_DB` to the `$VH2_ROOT` layout above. Defaults left on: three
assemblers, multi-source coverage, MultiQC on; QUAST / geNomad / dedup
off; `HOST_REMOVAL: bwa`.

## 5. Run across every run folder

`RESULT_FOLDER` = `RESULTS_FOLDER/<basename of SAMPLES>`, so each run
gets its own results subdirectory and its own `run_information_<run>.csv`
and `per_virus_<run>.csv`. Override `SAMPLES` per invocation.

```bash
conda activate virushanter
cd "$VH2_ROOT/code/virusHanter2"

# Validate discovery first (dry-run on one folder):
snakemake -n --sdm conda --configfile config/config.local.yaml \
  --config SAMPLES="$VH2_ROOT/runs/run01"

# Process every run folder. The first invocation also builds the per-rule
# conda envs (slower); later runs reuse them.
for run in "$VH2_ROOT"/runs/*/; do
  echo ">>> $run"
  snakemake --sdm conda --cores 16 \
    --configfile config/config.local.yaml \
    --config SAMPLES="$run" \
    || { echo "FAILED: $run"; break; }
done
```

Each run produces, under `$VH2_ROOT/results/<run>/`: per-sample HTML
reports, `run_information_<run>.csv`, `per_virus_<run>.csv`, mosdepth
summaries, and a batch `multiqc_report.html`. If a run is interrupted:
`snakemake --unlock --configfile config/config.local.yaml --config SAMPLES=<run>`
then re-run it.

## 6. Final merged output across all runs

`scripts/merge_runs.py` globs `run_information_*.csv` and
`per_virus_*.csv` from each result folder and concatenates them.

```bash
python scripts/merge_runs.py \
  $(for run in "$VH2_ROOT"/runs/*/; do printf ' --result-folder %s/results/%s' "$VH2_ROOT" "$(basename "$run")"; done) \
  --out-dir "$VH2_ROOT/master"
# writes:
#   $VH2_ROOT/master/master_per_sample.csv  (one row per sample; run_information schema)
#   $VH2_ROOT/master/master_per_virus.csv   (one row per sample x detected virus; the
#                                            16-column per-virus table)
```

`master_per_virus.csv` is the all-samples deliverable; the column
mapping to friendly collaborator labels is documented in
[PER_VIRUS_OUTPUT.md](PER_VIRUS_OUTPUT.md).

## Validation (optional)

To confirm the refactored pipeline reproduces the original
`virusHanter` for a batch, follow [PARITY_SIGNOFF.md](PARITY_SIGNOFF.md)
(use its parity-recovery config, e.g. `ASSEMBLERS: ["MEGAHIT"]`). For
production detection, leave the three-assembler defaults on.
