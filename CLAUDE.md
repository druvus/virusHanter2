# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working inside
`virusHanter2/`. The repository-wide CLAUDE.md at the regionen root
applies as well; this file adds the package-specific detail.

## Scope

`virusHanter2` is the data-processing half of the viral metagenomics
pipeline. It owns assembly, classification, alignment and
aggregation. HTML rendering is delegated to the
[`reportHanter`](../reportHanter) package and must not be
re-implemented here. The pairing must reproduce the parity-locked
report and `run_information_<batch>.csv` schema of the original
`virusHanter`; see [`docs/PARITY_NOTES.md`](docs/PARITY_NOTES.md).

## Layout

```
rules/
  pre_processing.smk    fastp, bwa host removal, markdup, optional dedup
  classification.smk    Kaiju, Kraken2, top-N viral hit selection
  assembly.smk          MEGAHIT + metaSPAdes (per-assembler wildcard),
                        Pilon, BLASTN, CheckV, optional QUAST and
                        geNomad
  post_processing.smk   mosdepth, generate_report, MultiQC, per-virus
                        and run-info aggregations
scripts/
  functions.py          Snakefile-side helpers (no plotting here)
  per_virus_metrics.py  Per-(sample, virus) CSV builder
  aggregate_*.py        Run-level concatenation
  build_virus_parquet.py  Reference parquet builder (utility)
config/                 Production, default and JSON schema
envs/                   One conda yaml per rule
refresh/                Standalone Snakemake workflow that rebuilds
                        VIRUS_PARQUET from NCBI Virus; not part of
                        the main pipeline DAG
test/                   Mini-fixture build_fixtures + run_smoke
docs/                   Long-form documentation
```

## Entry points

```bash
conda activate virushanter
snakemake -n --sdm conda --configfile config/config.local.yaml   # dry-run
snakemake    --sdm conda --cores N --configfile config/config.local.yaml
snakemake --unlock                                                # after Ctrl-C
```

The driver env only needs Snakemake (`snakemake-minimal=9.14.*` to
stay pickle-compatible with the `reporthanter` rule env). Every
per-rule tool gets its own conda env from `envs/*.yaml`, materialised
on first use.

## Optional stages — config-flag gated

| Flag | Default | Effect |
|---|---|---|
| `MULTIQC` | `"TRUE"` | Workflow-level MultiQC HTML at `{batch}/multiqc_report.html`. |
| `DEDUPLICATE` | `"FALSE"` | `remove_host` reads the markdup BAM with `-F 1024`; PCR duplicates are excluded from MEGAHIT and the coverage step. Off by default to preserve byte-identical parity. |
| `QUAST` | `"FALSE"` | `quast_per_assembler` runs against each (sample, assembler) pair and the per-assembler reports are fed to MultiQC. Bioconda has no `osx-arm64` build of QUAST; on Apple Silicon either keep this off or set `CONDA_SUBDIR=osx-64`. |
| `GENOMAD` | `"FALSE"` | `genomad` end-to-end runs alongside CheckV per (sample, assembler); geNomad's per-contig scores are appended as additive columns in `per_virus_metrics.csv`. Requires `GENOMAD_DB`. |
| `GENOMAD_SPLITS` | `4` | `genomad --splits N` value. Higher reduces peak `mmseqs prefilter` memory at the cost of run time; default 4 keeps the peak under ~6 GB on the DRRKK samples. Set `0` to restore mmseqs' auto-split. |
| `ASSEMBLERS` | `["MEGAHIT", "SPAdes"]` | List of de novo assemblers run per sample. Choices: `MEGAHIT`, `SPAdes` (metaSPAdes `--meta`), `rnaviralSPAdes` (SPAdes `--rnaviral`, tuned for RNA virus libraries). Each entry drives an independent Pilon / BLASTN / CheckV (and optional geNomad / QUAST) chain under `{sample}/{assembler}/...`. Defaults to MEGAHIT + metaSPAdes, which intentionally breaks byte-identical parity with the original `virusHanter` (see `docs/PARITY_NOTES.md`). Set `["MEGAHIT"]` to recover parity. |
| `HOST_REMOVAL` | `"bwa"` | Host-removal backend. `bwa` runs `bwa mem -k 26` against `HUMAN_INDEX` and is the parity default. `hostile` (Bede et al.) runs minimap2 against a bundled T2T-CHM13 + alt + decoy reference and catches telomeric, centromeric and segmental-duplication host reads that GRCh38-based BWA misses. `HOSTILE_INDEX` optionally points at a pre-downloaded T2T-CHM13 bundle. |
| `COVERAGE_SOURCES` | `["KRAKEN", "KAIJU", "BLAST"]` | Classifiers whose viral hits contribute taxids to the BWA reference set used by mosdepth coverage. The union of the top-N from each enabled source drives reference selection; an `unmapped_taxids.tsv` sidecar lists classified taxids missing from `VIRUS_PARQUET`. Set `["KRAKEN"]` to recover the pre-multi-source behaviour. |
| `COVERAGE_TOP_N` | `20` | Per-classifier cap on the number of viral hits whose taxids enter the BWA reference set. Applied independently to each entry in `COVERAGE_SOURCES`. |
| `TAXDUMP_NODES` | `""` | Optional path to an uncompressed NCBI `nodes.dmp`. When supplied, the coverage rule applies the rank filter and (when enabled) walks missing taxids up to their genus. Refresh alongside the parquet via `refresh/refresh_virus_parquet.smk`. |
| `COVERAGE_RANK_FILTER` | realm/kingdom/...family list | NCBI ranks that classifier hits are dropped at before they enter the union, so higher-rank propagation rows (Viruses, Cardeaviricetes, Anelloviridae, ...) no longer flood `unmapped_taxids.tsv`. Requires `TAXDUMP_NODES`; set to `[]` to disable. |
| `COVERAGE_GENUS_WALKUP` | `"TRUE"` | When a classifier hit is absent from `VIRUS_PARQUET`, walk up to its genus via the taxdump and substitute a genus reference. Tagged in the `virus_names` `sources` column with the `->genus` suffix. Requires `TAXDUMP_NODES`. |
| `COVERAGE_WINDOW` | `100` | Window size (bp) for `mosdepth --by`. Drives coverage-trace resolution in the HTML report. |

Every opt-in is parity-safe by default. When you flip one on, verify
the parity-locked columns of `run_information_<batch>.csv` are
unchanged; trailing additive columns are fine.

## Conventions

- Snakemake 9+, Python 3.12+ (matches the `panel` and `reporthanter`
  envs that share the driver tree).
- ASCII only in `*.smk` files; no Unicode in Snakemake or Nextflow
  files (project rule).
- Modest, plain scientific British English in documentation,
  docstrings and comments.
- Sample discovery is regex-based over `config["SAMPLES"]` — no
  sample sheet. `SECONDARY_HOST_OR_NOT` toggles an optional second
  host stage; several downstream rules use `lambda wildcards:`
  inputs that branch on it. Preserve that pattern when adding rules
  that consume host-removed reads.
- Trailing additions to `run_information_<batch>.csv` are fine;
  existing column positions and values must stay byte-identical to
  the original `virusHanter` output. Drop new trailing columns
  before diffing.
- Rules that invoke external tools use `script:` (or `shell:`),
  not `run:`. Snakemake's `conda:` directive is silently ignored
  for `run:` blocks, so a `run:` rule that calls a binary from
  inside Python ends up needing the binary in the *driver* env
  rather than the per-rule env. The tool-running rules live in
  `scripts/run_<rule>.py` and take the standard `snakemake`
  magic object. Pandas-only rules (e.g. `wrangle_pilon`,
  `merge_checkv_blastn`) may stay as `run:` because pandas is in
  the driver env.
- Every classifier and BLAST output is canonicalised to the
  ICTV-binomial species name via the NCBI taxdump's
  `find_species_taxid` walk-up; the rewriter also appends an
  `aliases` column carrying the legacy NCBI scientific name plus
  the `names.dmp` alias categories (acronym, common name,
  equivalent name, ...). See
  [docs/PIPELINE.md#canonical-species-naming-ictv-binomials](docs/PIPELINE.md#canonical-species-naming-ictv-binomials).
  When adding a new rule that emits a per-tax_id column, run it
  through `canonicalise_taxon_names` so the report stays
  consistent.

## Testing

`test/run_smoke.sh` runs lint + DAG dry-run by default; `--full`
builds tiny fixtures under `test/mini_db/` and runs the pipeline
through `generate_report` if a CheckV database is present. The
fixture has three samples (`sample1_R`, `sample2_R`, `sample3_R`),
each carrying a distinct synthetic virus (`alpha` / `beta` /
`gamma`, taxids 100001-100003), so the smoke exercises multi-sample
paths and the per-virus attribution. See
[`test/README.md`](test/README.md) for the details.

Apple Silicon caveats:

- MEGAHIT's `_no_hw_accel count -k 21` SIGSEGVs on small inputs;
  `rules/assembly.smk` already passes `--k-min 27` on
  `Darwin/arm64` to dodge the buggy path. The retry loop (default
  4 attempts via `MEGAHIT_RETRIES`) handles the residual
  non-determinism.
- A CheckV DB sourced from an external Mac volume can carry
  AppleDouble `._*.hmm` metadata files that CheckV mis-reads as
  HMMs ("80 hmmsearch tasks failed"). Strip them with
  `find /path/to/checkv-db-vX -name '._*' -delete` once.
- geNomad's embedded `mmseqs prefilter` can OOM on
  metaSPAdes-derived proteomes (3× larger than MEGAHIT). The
  default `GENOMAD_SPLITS: 4` partitions the search; raise it on
  smaller hosts.
- `wget` on a LaCie-mounted external volume emits a non-fatal
  `utime()` warning that some builds expose as a non-zero exit
  code. The refresh workflow uses `curl` for the large downloads
  to avoid the issue; mirror that pattern in any new download
  rule.

When the CheckV DB is stubbed, the smoke degrades automatically to
`--until blastn mosdepth_kraken_hits kaiju_to_table`.

## Refresh workflow

`refresh/refresh_virus_parquet.smk` is a standalone Snakemake
workflow (not part of the main pipeline DAG) that rebuilds
`VIRUS_PARQUET`, the matching Kaiju FMI index, the Kraken2 viral
DB and the BLAST viral DB tarballs from the same NCBI viral
RefSeq snapshot, downloads the taxdump, and emits an
overlap-with-Kraken2 diagnostic sidecar. Building Kraken2 from
the same RefSeq pull closes the recurring gap where the publicly
hosted `k2_viral_*` snapshots occasionally omit individual
genomes (e.g. the Feb 2026 snapshot missed HSV-2 / NC_001798).
Driven by `refresh/config.yaml`. Helpers live at
`scripts/build_virus_parquet.py`,
`scripts/reformat_kaiju_headers.py`,
`scripts/compare_parquet_kraken2.py`. See
[`docs/REFRESH_TUTORIAL.md`](docs/REFRESH_TUTORIAL.md).
