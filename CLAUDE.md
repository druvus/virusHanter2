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
  assembly.smk          MEGAHIT, Pilon, BLASTN, CheckV, optional QUAST
                        and geNomad
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
| `QUAST` | `"FALSE"` | `quast_megahit` runs against each sample's MEGAHIT contigs and is fed to MultiQC. Bioconda has no `osx-arm64` build of QUAST; on Apple Silicon either keep this off or set `CONDA_SUBDIR=osx-64`. |
| `GENOMAD` | `"FALSE"` | `genomad` end-to-end runs alongside CheckV; geNomad's per-contig scores are appended as additive columns in `per_virus_metrics.csv`. Requires `GENOMAD_DB`. |
| `ASSEMBLERS` | `["MEGAHIT", "SPAdes"]` | List of de novo assemblers run per sample. Each entry drives an independent Pilon / BLASTN / CheckV (and optional geNomad / QUAST) chain under `{sample}/{assembler}/...`. Defaults to both assemblers, which intentionally breaks byte-identical parity with the original `virusHanter` (see `docs/PARITY_NOTES.md`). Set `["MEGAHIT"]` to recover parity. |
| `COVERAGE_SOURCES` | `["KRAKEN", "KAIJU", "BLAST"]` | Classifiers whose viral hits contribute taxids to the BWA reference set used by mosdepth coverage. The union of the top-N from each enabled source drives reference selection; an `unmapped_taxids.tsv` sidecar lists classified taxids missing from `VIRUS_PARQUET`. Set `["KRAKEN"]` to recover the pre-multi-source behaviour. |
| `COVERAGE_TOP_N` | `20` | Per-classifier cap on the number of viral hits whose taxids enter the BWA reference set. Applied independently to each entry in `COVERAGE_SOURCES`. |
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
  `Darwin/arm64` to dodge the buggy path.
- A CheckV DB sourced from an external Mac volume can carry
  AppleDouble `._*.hmm` metadata files that CheckV mis-reads as
  HMMs ("80 hmmsearch tasks failed"). Strip them with
  `find /path/to/checkv-db-vX -name '._*' -delete` once.

When the CheckV DB is stubbed, the smoke degrades automatically to
`--until blastn mosdepth_kraken_hits kaiju_to_table`.
