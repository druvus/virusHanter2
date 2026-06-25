# Configuration

The pipeline reads two files at startup: `config/config.yaml` (your
inputs and database paths) and `config/config.schema.yaml` (a JSON
schema that fails fast on placeholder paths). All Conda environments
are declared in `envs/*.yaml` and materialised on demand by
Snakemake.

## `config/config.yaml`

### Required keys

| Key | Purpose |
|---|---|
| `SAMPLES` | Directory containing paired-end FASTQs. Sample names are derived from the common prefix of each `_R1`/`_R2` pair; both `.fastq` and `.fastq.gz` inputs are picked up. |
| `RESULTS_FOLDER` | Output root. The actual results live at `{RESULTS_FOLDER}/{basename(SAMPLES)}/`. |
| `THREADS` | Default per-rule thread count; combined with `snakemake --cores N` for scheduling. |
| `HUMAN_INDEX` | BWA index prefix for the host (human) genome. |
| `KAIJU_DB` | Directory containing `.fmi`, `names.dmp`, `nodes.dmp`. |
| `KRAKEN_DB` | Kraken2 database directory (with `hash.k2d`, `opts.k2d`, `taxo.k2d`). |
| `BLASTN_DB` | BLAST nucleotide database **prefix** (e.g. the `.nal` alias name). |
| `CHECKV_DB` | CheckV database directory. |
| `VIRUS_PARQUET` | Parquet with columns `(name, sequence, tax_id, rank, genus_taxid)` used by `bwa_align_to_kraken_hits`. The `rank` and `genus_taxid` columns are added when the parquet is built with `--taxdump-nodes`; older 3-column parquets still work but the rank filter + genus walk-up degrade to no-ops. See [REFERENCE_DBS.md](REFERENCE_DBS.md) and the [refresh tutorial](REFRESH_TUTORIAL.md). |
| `NUMBER_OF_PLOTS` | Top-N Kraken viral hits surfaced in the per-virus CSV (default `10`). |

### Always-on stages, with tuning knobs

| Key | Default | Purpose |
|---|---|---|
| `CONTIG_LENGTH` | `500` | Minimum polished-contig length kept after Pilon. |
| `PILON_MEM` | `"50G"` | JVM heap for Pilon. |
| `MEGAHIT_MEM_FRACTION` | `0.5` | Fraction of system RAM MEGAHIT is allowed to allocate. Drop to `0.2-0.3` on memory-tight laptops. |
| `MEGAHIT_RETRIES` | `4` | Apple Silicon MEGAHIT non-determinism mitigation: number of retry attempts before falling back to the `DUMMY_CONTIG`. Linux defaults to a single attempt. |
| `COVERAGE_WINDOW` | `100` | Window size (bp) passed to `mosdepth --by`. Smaller values give a finer coverage trace in the report and a larger `regions.bed.gz`. |

### Optional stages (config-flag gated)

The table below mirrors the opt-in stage table in
[`CLAUDE.md`](../CLAUDE.md). Every flag is parity-safe at its
default; flipping it to its non-default value is a deliberate
divergence from the original `virusHanter`.

| Flag | Default | Effect |
|---|---|---|
| `CLEAN` | `"FALSE"` | If `"TRUE"`, remove intermediates after the run and write `analysis_done.txt`. |
| `MULTIQC` | `"TRUE"` | Emit `{batch}/multiqc_report.html` at the end of the run. |
| `DEDUPLICATE` | `"FALSE"` | Exclude PCR duplicates from the host-removed reads that feed assembly and the coverage step. |
| `QUAST` | `"FALSE"` | Run QUAST per (sample, assembler) and feed each report to MultiQC. Bioconda has no `osx-arm64` build of QUAST. |
| `GENOMAD` | `"FALSE"` | Run geNomad alongside CheckV per (sample, assembler). Requires `GENOMAD_DB`. |
| `GENOMAD_DB` | `""` | Path to the populated `genomad_db/` directory. |
| `GENOMAD_SPLITS` | `4` | `genomad --splits N` value. Higher reduces peak `mmseqs prefilter` memory at the cost of run time; default 4 keeps the peak under ~6 GB. Set `0` on hosts with abundant RAM to restore mmseqs' auto-split. |
| `ASSEMBLERS` | `["MEGAHIT", "metaSPAdes", "rnaviralSPAdes"]` | List of de novo assemblers run per sample. Choices: `MEGAHIT`, `metaSPAdes` (SPAdes `--meta`), `rnaviralSPAdes` (SPAdes `--rnaviral`, tuned for RNA virus libraries). Defaults to all three; the contigs of all three fold into `number_of_contigs` and `top_contigs_blastn`. Each entry drives an independent Pilon / BLASTN / CheckV (and optional geNomad / QUAST) chain under `{sample}/{assembler}/`. Set to `["MEGAHIT"]` for parity with the original `virusHanter`. The deprecated alias `SPAdes` is rejected at workflow load; use `metaSPAdes`. |
| `HOST_REMOVAL` | `"bwa"` | Host-removal backend. `bwa` runs `bwa mem -k 26` against `HUMAN_INDEX` (parity default). `hostile` (Bede et al.) runs minimap2 against a bundled T2T-CHM13 + alt + decoy human reference and catches telomeric, centromeric and segmental-duplication host reads that GRCh38-based BWA misses. The active backend's name lands in the trailing `host_removal_tool` column of `run_information_<batch>.csv`. |
| `HOSTILE_INDEX` | `""` | Optional override for the hostile index (path to a pre-downloaded directory or a bare hostile-managed name). Only consulted when `HOST_REMOVAL: "hostile"`. When empty, the rule passes `--index human-t2t-hla.rs-viral-202401_ml-phage-202401` — the masked variant that keeps endogenous retrovirus, phage and other host-embedded viral reads (the bare `human-t2t-hla` index would drop them as host). |
| `COVERAGE_SOURCES` | `["KRAKEN", "KAIJU", "BLAST"]` | Classifiers whose viral hits contribute tax_ids to the BWA reference set used by mosdepth coverage. Union of the per-classifier top-N drives reference selection. Set to `["KRAKEN"]` to recover the pre-multi-source behaviour. |
| `COVERAGE_TOP_N` | `20` | Per-classifier cap on the number of viral hits whose tax_ids enter the BWA reference set. |
| `TAXDUMP_NODES` | `""` | Optional path to an uncompressed NCBI `nodes.dmp`. When set, enables the **ICTV species walk-up** that rewrites every classifier and BLAST output to the species-rank scientific name and adds an `aliases` column with the legacy NCBI names (see [Canonical species naming](PIPELINE.md#canonical-species-naming-ictv-binomials)); plus the rank filter and (when on) the genus walk-up. The sibling `names.dmp` must live next to `nodes.dmp` (the refresh workflow publishes both). Built and published as part of the refresh workflow; see [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md). |
| `COVERAGE_RANK_FILTER` | `[acellular root, realm, kingdom, subkingdom, phylum, subphylum, class, subclass, order, suborder, family, subfamily]` | NCBI rank strings that classifier hits are dropped at before they enter the coverage union. Higher-rank propagation rows have no per-tax_id sequence so they would otherwise flood `unmapped_taxids.tsv`. Requires `TAXDUMP_NODES`; set to `[]` to disable. |
| `COVERAGE_GENUS_WALKUP` | `"TRUE"` | When a classifier hit is absent from `VIRUS_PARQUET`, walk up to its genus and substitute a representative genus reference. Tagged in the `virus_names` `sources` column with the `->genus` suffix. Requires `TAXDUMP_NODES`. |
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
| `envs/samtools.yaml` | samtools |
| `envs/kraken.yaml` | kraken2 |
| `envs/kaiju.yaml` | kaiju |
| `envs/megahit.yaml` | megahit |
| `envs/spades.yaml` | spades / metaspades |
| `envs/pilon.yaml` | pilon, bwa, samtools, openjdk |
| `envs/blastn.yaml` | blast, pandas, pyfastx |
| `envs/checkv.yaml` | checkv |
| `envs/mosdepth.yaml` | mosdepth |
| `envs/multiqc.yaml` | multiqc |
| `envs/quast.yaml` | quast (only used when `QUAST: "TRUE"`) |
| `envs/genomad.yaml` | genomad |
| `envs/panel.yaml` | python, pandas, pyfastx, pyarrow (wrangling rules) |
| `envs/reporthanter.yaml` | python 3.12 + pip-installed `reportHanter` from GitHub |
| `envs/refresh.yaml` | python + pandas + pyarrow + pyfastx + curl + wget + kraken2 + kaiju (used only by the refresh workflow under `refresh/`) |
| `envs/hostile.yaml` | hostile + minimap2 + samtools (only used when `HOST_REMOVAL: "hostile"`) |

Snakemake materialises an env the first time any rule that declares
it runs. Subsequent runs reuse the cached env under
`.snakemake/conda/<hash>/`.

## Apple Silicon notes

The driver env stays native `osx-arm64`. A few bioconda tools have
historical or current rough edges on this platform; see
[`../test/run_smoke.sh`](../test/run_smoke.sh) for the explicit
workarounds. The current pipeline already encodes them:

- MEGAHIT runs single-threaded with `--no-hw-accel --k-min 27
  --k-max 57`, and the retry loop catches the residual SIGSEGV /
  SIGABRT non-determinism.
- `GENOMAD_SPLITS: 4` keeps the embedded `mmseqs prefilter` under
  the 18 GB system RAM ceiling typical of an Apple Silicon laptop.
- QUAST has no `osx-arm64` bioconda build; either leave
  `QUAST: "FALSE"` or set `CONDA_SUBDIR=osx-64`.
