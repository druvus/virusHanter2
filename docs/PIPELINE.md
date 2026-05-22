# Pipeline

The workflow is described as Snakemake rules under `rules/`,
included into a single `Snakefile`. Rules are grouped by stage.
Throughout, `{assembler}` is a wildcard that takes a value from
`config[ASSEMBLERS]` (default `["MEGAHIT", "SPAdes"]`), so every
contig-producing rule below the de novo assemblers runs once per
(sample, assembler) pair and lands its outputs under
`{batch}/{sample}/{assembler}/...`.

## Stages

### Pre-processing — `rules/pre_processing.smk`

1. **`fastp`** — quality trim/filter the paired-end input.
2. **Host removal** — one of two backends, selected by
   `config[HOST_REMOVAL]`:
   - **`bwa`** (default, parity-safe): `bwa_human` →
     `markdup_human` (markdup stats only) → `remove_host`
     (flagstat + `samtools view -f 12` to extract pairs where
     both reads are unmapped to the human host) →
     `bam_to_fastq_human` (paired FASTQ).
   - **`hostile`**: `hostile_human` runs Bede et al.'s hostile
     against the bundled T2T-CHM13 reference (more thorough
     for telomeric / pericentromeric host reads). Produces the
     same paired FASTQ outputs plus a flagstat-shape stats
     file the report layer parses.
3. Downstream rules consume the host-removed FASTQ via
   `host_removed_r1` / `host_removed_r2` helper functions that
   resolve to the active backend's outputs.
4. **`bwa_secondary_host`** + **`remove_secondary_host`** +
   **`bam_to_fastq_secondary`** — optional second host removal
   when `SECONDARY_HOST_INDEX` is set.

### Classification — `rules/classification.smk`

7. **`kaiju`** + **`kaiju_to_table`** — protein-level taxonomic
   classification against the configured `KAIJU_DB`.
8. **`kraken`** + **`wrangle_kraken`** — *k*-mer DNA classification
   against the configured `KRAKEN_DB`; the wrangled CSV adds an
   explicit `domain` column.

### Assembly + annotation — `rules/assembly.smk`

Every step below runs once per assembler in `config[ASSEMBLERS]`.
Outputs land under `{sample}/{assembler}/...`.

9. **`megahit`** — de novo assembly. Apple Silicon retry loop
   (`MEGAHIT_RETRIES`) handles the bioconda `osx-arm64`
   non-determinism; falls back to a `DUMMY_CONTIG` when every
   retry crashes so downstream rules always have an input.
10. **`metaspades`** — de novo assembly via metaSPAdes
    (`--meta`). Mirrors the dummy-contig fallback on assembly
    refusal.
11. **`rnaviralspades`** *(optional)* — SPAdes `--rnaviral` for
    RNA-virus-tuned assembly. Activated by adding
    `"rnaviralSPAdes"` to `config[ASSEMBLERS]`.
11. **`pilon`** + **`wrangle_pilon`** — short-read polishing, then
    length-filter to `CONTIG_LENGTH`. `wrangle_pilon` stamps an
    `assembler` column on the per-contig CSV so every downstream
    consumer knows which assembler produced each row.
12. **`blastn`** — best-hit annotation against `BLASTN_DB`.
13. **`checkv`** — viral contig contamination / completeness call.
14. **`merge_checkv_blastn`** — inner join CheckV columns into the
    BLASTN table. The merged CSV (one per (sample, assembler)) is
    what the per-sample HTML report and the per-virus aggregator
    consume.
15. **`genomad`** *(optional)* — second viral-contig classifier
    when `GENOMAD: "TRUE"`. One geNomad run per (sample, assembler);
    the `--splits` flag is configurable via `GENOMAD_SPLITS`
    (default 4) to keep peak `mmseqs prefilter` memory bounded.
16. **`quast_per_assembler`** *(optional)* — N50, largest contig,
    GC%, total length per (sample, assembler) when `QUAST: "TRUE"`.

### Post-processing + reporting — `rules/post_processing.smk`

17. **`bwa_align_to_kraken_hits`** — build a per-sample reference
    set from the **union** of every classifier in
    `config[COVERAGE_SOURCES]`:
    - Kraken: top-`COVERAGE_TOP_N` viral taxa.
    - Kaiju: top-`COVERAGE_TOP_N` viral taxa filtered against the
      parquet's tax_id set.
    - BLAST: every per-assembler merged CSV; accessions are
      resolved to tax_ids via the parquet.
    When `TAXDUMP_NODES` is configured, a **rank filter**
    (`COVERAGE_RANK_FILTER`) drops higher-rank propagation rows
    (kingdom / phylum / class / order / family) silently, and a
    **genus walk-up** (`COVERAGE_GENUS_WALKUP`) substitutes a
    representative genus reference when the exact tax_id is absent
    from the parquet. Emits the BWA index, the per-reference
    `kraken_top_viruses.fasta`, a `virus_names` sidecar with a
    `sources` column, and an `unmapped_taxids.tsv` audit trail
    listing classifier hits that could not be served.
18. **`mosdepth_kraken_hits`** — per-reference coverage stats with
    `--by COVERAGE_WINDOW` and `--thresholds 1,5,10`. The
    `regions.bed.gz` is consumed by `reporthanter` to render
    interactive coverage traces per reference.
19. **`generate_report`** — invoke `reporthanter` v0.5+ to render
    the per-sample interactive HTML. Inputs include one BLAST
    merged CSV per active assembler (and one QUAST report / one
    geNomad summary per assembler when those flags are on).
20. **`aggregate_run_information`** — concatenate per-sample
    summaries into `run_information_<batch>.csv`. Trailing
    per-assembler columns (`<assembler>_n_contigs`,
    `<assembler>_n50`, plus `assemblers_used`) sit beside the
    parity-locked legacy columns.
21. **`per_virus_metrics`** + **`aggregate_per_virus`** — join
    Kraken / Kaiju / BLASTN / mosdepth into per-(sample, virus)
    rows and write the collaborator-facing `per_virus_<batch>.csv`.
    The per-assembler `<assembler>_contigs` counts are trailing
    additive columns; schema in [PER_VIRUS_OUTPUT.md](PER_VIRUS_OUTPUT.md).
22. **`multiqc`** — workflow-level QC dashboard (gated by
    `MULTIQC: "TRUE"`).
23. **`clean_everything`** — optional cleanup of intermediates when
    `CLEAN: "TRUE"`.

## Outputs

Everything lives under `{RESULTS_FOLDER}/{batch}/`, where `batch`
is the basename of `SAMPLES`.

### Per-sample (`{batch}/{sample}/`)

| Path | Contents |
|---|---|
| `FASTP/<sample>.fastp.{html,json}` | fastp QC report. |
| `FASTP/<sample>_r{1,2}_trimmed.fq` | trimmed reads. |
| `bwa/<sample>_human.bam`, `_human_unmapped.bam`, `_human_unmapped_r{1,2}.fastq` | host alignment + unmapped read extraction. |
| `logs/human_contamination_flagstat.txt` | flagstat over the host alignment. |
| `logs/human_markdup_stats.txt` | `samtools markdup -s` summary. |
| `KAIJU/<sample>.kaiju.{out,table.tsv}` | Kaiju classification. |
| `KRAKEN/<sample>.kraken.{report,csv}` | Kraken2 report + wrangled CSV. |
| `{ASSEMBLER}/<sample>.contigs.fa` | de novo contigs, one folder per active assembler. |
| `{ASSEMBLER}/PILON/<sample>_improved_contigs.fasta`, `<sample>.contigs.csv` | polished + length-filtered. |
| `{ASSEMBLER}/BLASTN/<sample>.contigs.blastn.csv` | best-hit annotation. |
| `{ASSEMBLER}/CHECKV/<sample>.contamination.tsv`, `<sample>.merged.csv` | CheckV call + BLASTN+CheckV inner join. |
| `{ASSEMBLER}/GENOMAD/<sample>_improved_contigs_summary/<sample>_improved_contigs_virus_summary.tsv` | only when `GENOMAD: "TRUE"`. |
| `{ASSEMBLER}/QUAST/report.tsv` | only when `QUAST: "TRUE"`. Rendered as a sub-tab per assembler under the report's dedicated "Assembly" section (sitting between Classification of Raw Reads and Classification of Contigs in the data-flow order assembly → annotation), since QUAST measures the assembler's contigs rather than the host alignment. |
| `BWA_KRAKEN/<sample>_kraken.bam`, `kraken_top_viruses.fasta`, `kraken_top_virus_names.tsv`, `unmapped_taxids.tsv` | multi-source reference set, sidecar with `sources` column, audit list of classifier hits without a parquet reference. |
| `MOSDEPTH/<sample>.mosdepth.summary.txt`, `.regions.bed.gz`, `.thresholds.bed.gz` | numeric coverage (drives the per-reference traces in the HTML report). |
| `REPORT/<sample>.html` | per-sample interactive report. |
| `<sample>.per_virus.csv` | one row per detected Kraken viral taxid, with trailing per-assembler contig counts. |

### Per-batch (`{batch}/`)

| Path | Contents |
|---|---|
| `run_information_<batch>.csv` | one row per sample; parity-locked legacy columns + trailing per-assembler stats. |
| `per_virus_<batch>.csv` | concatenation of every per-sample per-virus CSV. |
| `multiqc_report.html` + `multiqc_data/` | run-level QC dashboard. |
| `analysis_done.txt` | sentinel when `CLEAN: "TRUE"`. |

### Combining multiple Illumina runs

`scripts/merge_runs.py` is a standalone CLI: pass a `--result-folder`
per batch, get back `master_per_sample.csv` and `master_per_virus.csv`.
See [PER_VIRUS_OUTPUT.md](PER_VIRUS_OUTPUT.md#multi-run-master-files)
for the full schema and a worked example.
