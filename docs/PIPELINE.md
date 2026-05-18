# Pipeline

The workflow is described as Snakemake rules under `rules/`,
included into a single `Snakefile`. Rules are grouped by stage.

## Stages

### Pre-processing — `rules/pre_processing.smk`

1. **`fastp`** — quality trim/filter the paired-end input.
2. **`bwa_human`** — map trimmed reads to the human host (`bwa mem -k 26`).
3. **`remove_host`** — `samtools flagstat` + extract unmapped pairs.
4. **`bam_to_fastq_human`** — write host-unmapped pairs back to FASTQ.
5. **`markdup_human`** — `samtools markdup -s` on the host BAM. Stats-only side channel; does not change the un-marked BAM that downstream rules read.
6. **`bwa_secondary_host`** + **`remove_secondary_host`** + **`bam_to_fastq_secondary`** — optional second host removal when `SECONDARY_HOST_INDEX` is set.

### Classification — `rules/classification.smk`

7. **`kaiju`** + **`kaiju_to_table`** — protein-level taxonomic classification.
8. **`kraken`** + **`wrangle_kraken`** — k-mer DNA classification; the wrangled CSV adds an explicit `domain` column.

### Assembly + annotation — `rules/assembly.smk`

9. **`megahit`** — de novo assembly. Falls back to a `DUMMY_CONTIG` when the assembler emits nothing, so downstream rules always have an input.
10. **`pilon`** + **`wrangle_pilon`** — short-read polishing, then length-filter to `CONTIG_LENGTH` and emit a CSV of polished contigs.
11. **`blastn`** — best-hit annotation against the configured viral nucleotide DB.
12. **`checkv`** — viral contig contamination / completeness call.
13. **`merge_checkv_blastn`** — inner join CheckV columns into the BLASTN table; this is what the per-sample HTML report consumes.
14. **`genomad`** *(optional)* — second viral-contig classifier when `GENOMAD: "TRUE"`. Writes a per-sample summary TSV under `GENOMAD/`. Does not feed the report; sits alongside CheckV.

### Post-processing + reporting — `rules/post_processing.smk`

15. **`bwa_align_to_kraken_hits`** — pick the top-N Kraken viral taxa, look up their reference sequences in `VIRUS_PARQUET`, build a BWA index of just those references, and map the host-unmapped reads.
16. **`bam2plot`** — SVG coverage profiles for references that exceed `PLOT_THRESHOLD`.
17. **`mosdepth_kraken_hits`** — numeric per-reference coverage stats with `--thresholds 1,5,10`.
18. **`generate_report`** — invoke `reporthanter` to render the per-sample interactive HTML.
19. **`aggregate_run_information`** — concatenate per-sample summaries into `run_information_<batch>.csv` (parity-locked to the original `virusHanter`; see [PARITY_NOTES.md](PARITY_NOTES.md)).
20. **`per_virus_metrics`** + **`aggregate_per_virus`** — join Kraken/Kaiju/BLASTN/mosdepth into per-(sample, virus) rows and write the collaborator-facing `per_virus_<batch>.csv`. Schema: [PER_VIRUS_OUTPUT.md](PER_VIRUS_OUTPUT.md).
21. **`multiqc`** — workflow-level QC dashboard (gated by `MULTIQC: "TRUE"`).
22. **`clean_everything`** — optional cleanup of intermediates when `CLEAN: "TRUE"`.

## Outputs

Everything lives under `{RESULTS_FOLDER}/{batch}/`, where `batch` is
the basename of `SAMPLES`.

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
| `MEGAHIT/<sample>.contigs.fa` | de novo contigs. |
| `PILON/<sample>_improved_contigs.fasta`, `<sample>.contigs.csv` | polished + length-filtered. |
| `BLASTN/<sample>.contigs.blastn.csv` | best-hit annotation. |
| `CHECKV/<sample>.contamination.tsv`, `<sample>.merged.csv` | CheckV call + BLASTN+CheckV inner join. |
| `GENOMAD/<sample>_summary/<sample>_virus_summary.tsv` | only when `GENOMAD: "TRUE"`. |
| `BWA_KRAKEN/<sample>_kraken.bam` + `kraken_top_viruses.fasta` | mapping to top-N Kraken viral references. |
| `COVERAGE_PLOTS/*.svg` | per-reference coverage profiles. |
| `MOSDEPTH/<sample>.mosdepth.summary.txt`, `.regions.bed.gz`, `.thresholds.bed.gz` | numeric coverage. |
| `REPORT/<sample>.html` | per-sample interactive report. |
| `<sample>.per_virus.csv` | one row per detected Kraken viral taxid. |

### Per-batch (`{batch}/`)

| Path | Contents |
|---|---|
| `run_information_<batch>.csv` | one row per sample, parity-locked schema. |
| `per_virus_<batch>.csv` | concatenation of every per-sample per-virus CSV. |
| `multiqc_report.html` + `multiqc_data/` | run-level QC dashboard. |
| `analysis_done.txt` | sentinel when `CLEAN: "TRUE"`. |

### Combining multiple Illumina runs

`scripts/merge_runs.py` is a standalone CLI: pass a `--result-folder`
per batch, get back `master_per_sample.csv` and `master_per_virus.csv`.
See [PER_VIRUS_OUTPUT.md](PER_VIRUS_OUTPUT.md#multi-run-master-files)
for the full schema and a worked example.
