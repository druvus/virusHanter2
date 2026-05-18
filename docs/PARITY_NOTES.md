# Parity notes: virusHanter -> virusHanter2 + reportHanter

The refactor splits the original monolithic `virusHanter/Snakefile` into a
data-processing pipeline (`virusHanter2`) and a report-generation package
(`reportHanter`). The two together must produce the same per-sample HTML
report and the same `run_information_<batch>.csv` schema as `virusHanter`.
This document captures the procedure used to validate that, and any
intended divergences that are expected to show up in a diff.

## Procedure

Pick a representative real sample (or, ideally, a batch of at least two
samples) that has been processed by the original `virusHanter` recently
enough that all required databases on the workstation are still consistent.

1. Run the original pipeline against the chosen samples folder using
   `virusHanter/Snakefile` and the database paths in
   `virusHanter/config.yaml`. Record the resulting `RESULTS_FOLDER`.

2. Run the refactored pipeline against the same samples folder using
   `virusHanter2/Snakefile` and a `config.yaml` that points at the same
   databases. Direct its output to a separate `RESULTS_FOLDER`.

   ```
   conda activate virushanter
   snakemake --use-conda --cores <N> --configfile config/config.yaml
   ```

3. Compare outputs for each sample.

### Per-sample HTML

The HTML reports are not byte-comparable (Panel embeds build-time
timestamps and randomized DOM ids). Compare instead by opening both files
side-by-side and checking each tab:

- Alignment Stats: total reads, percent mapped to human (and to the
  secondary host if configured), fastp summary rows.
- Classification of Raw Reads: top viral hits in the Kraken (virus and
  all-domains) panes and in the Kaiju pane; the unclassified percentages
  shown in the axis title.
- Classification of Contigs: contig table rows (name, match, accession,
  percent identity, sequence length).
- Alignment Coverage: same SVGs present.

### Run-aggregation CSV

`run_information_<batch>.csv` is byte-stable apart from the
`html_report` column (hex-encoded HTML, which carries the same
timestamp/build-id non-determinism described above). The
`kaiju_report` and `blastn_report` blob columns are byte-stable across
runs (no embedded timestamps). Diff the two CSVs with the HTML blob
column dropped, e.g.:

```
csvcut -C html_report run_information_<batch>.csv | sort > a.csv
csvcut -C html_report run_information_<batch>.csv | sort > b.csv
diff a.csv b.csv
```

Numeric columns (read counts, percentages, `number_of_contigs`) should
match exactly. String columns (`top_contigs_blastn`, `top_virus_kaiju`)
should match exactly unless one of the documented divergences below
applies. Column order should also match: `run_name, sample_name, date,
read_len, number_reads, mapped_to_human_percent, kraken_virus_percent,
kaiju_virus_percent, number_of_contigs, top_contigs_blastn,
top_virus_kaiju, html_report, kaiju_report, blastn_report`.

Two trailing columns (`duplicate_pairs`, `duplicate_rate_percent`)
were added on 2026-05-17 from the new `markdup_human` rule. They are
blank for any sample whose folder does not contain
`logs/human_markdup_stats.txt`, so dropping them from the new CSV
before diffing against an older one keeps the diff clean. See the
"Additive 2026-05-17 audit changes" section below.

## Recent parity work

The 2026-05-16 audit pass against `virusHanter/Snakefile` corrected
several silent numerical and string divergences that would have
surfaced as differences in a real-sample diff. Headline items:

- `aggregate_run_information.py` now appends `(reads)` / `(read_len)`
  suffixes to the `top_virus_kaiju` and `top_contigs_blastn` strings,
  applies `.dropna()` before summing the kaiju percent, restores the
  `run_name`, `date`, `kaiju_report`, and `blastn_report` columns, and
  uses the original column names (`read_len`, `html_report`).
- `bwa_human` re-adds `-k 26` to the bwa-mem invocation.
- `merge_checkv_blastn` reverts to an inner join on `name`.
- `bwa_align_to_kraken_hits` sorts viral taxa by `percent` descending
  before selecting the top 20.
- `scripts/functions.py:run_blastn` sets `BLASTDB` to the database
  parent directory so blastn can find any `taxdb.*` auxiliary files.
- `config/config.yaml` defaults align with the original: `CONTIG_LENGTH`
  is 500 and `PLOT_THRESHOLD` is 5.
- `reportHanter` Kraken `filter_data` default cutoff is 0.001 (matching
  the original `plot_kraken`).

## Additive 2026-05-17 audit changes

Following a targeted audit for the Twist Comprehensive Virus Research
Panel use case, three workflow-level additions were made. All are
purely additive: existing rules, per-sample tabs, and the original
fourteen `run_information_<batch>.csv` columns remain byte-identical
to a pre-change run.

- **`markdup_human` rule** (`rules/pre_processing.smk`): runs
  `samtools sort -n` -> `samtools fixmate -m` -> `samtools sort` ->
  `samtools markdup -s` on the bwa_human BAM and writes
  `logs/human_markdup_stats.txt`. The original bwa_human BAM is not
  modified; `remove_host` and every downstream rule continue to read
  the un-marked file, so flagstat-derived columns are unchanged.
- **`mosdepth_kraken_hits` rule** (`rules/post_processing.smk`):
  produces per-reference coverage summaries
  (`MOSDEPTH/<sample>.mosdepth.summary.txt` and
  `.regions.bed.gz`) from the same BAM `bwa_align_to_kraken_hits`
  emits. Sits alongside `bam2plot`'s `COVERAGE_PLOTS/` SVGs in a
  separate directory; bam2plot's outputs are unchanged.
- **`multiqc` rule** (`rules/post_processing.smk`): workflow-level
  rule that emits `{RESULT_FOLDER}/multiqc_report.html`. Gated by the
  `MULTIQC: "TRUE"` config flag (default on); set to `"FALSE"` to skip.
- **`per_virus_metrics` + `aggregate_per_virus` rules**
  (`rules/post_processing.smk`): join Kraken/Kaiju/BLASTN/mosdepth/
  fastp/flagstat into a flat per-(sample, virus) CSV. Per-sample
  output at `{sample}/{sample}.per_virus.csv`; batch-level
  concatenation at `{RESULT_FOLDER}/per_virus_<batch>.csv`. Schema
  documented in `docs/PER_VIRUS_OUTPUT.md`. Does not touch the
  existing per-sample summary CSV or HTML reports.
- **`mosdepth_kraken_hits` rule**: extended with `--thresholds 1,5,10`
  to emit `MOSDEPTH/<sample>.thresholds.bed.gz`. Pure addition;
  `summary.txt` and `regions.bed.gz` columns are unchanged from the
  original mosdepth invocation when callers ignore the threshold
  data.
- **`genomad` rule** (`rules/assembly.smk`): optional second
  viral-contig classifier. Off by default
  (`GENOMAD: "FALSE"`); when enabled, runs `genomad end-to-end` on
  the Pilon-polished contigs and writes a per-sample summary TSV
  under `GENOMAD/<sample>_summary/`. Does not merge into
  `merged_csv` or feed the report, so existing tabs/columns stay
  byte-identical regardless of the flag. Practical note: the
  bundled mmseqs2 step is heavily I/O-bound, so the geNomad DB
  must live on local SSD or fast networked storage — a USB-attached
  DB stalls for hours. See `docs/REFERENCE_DBS.md`.

Aggregated CSV: the two new columns `duplicate_pairs` and
`duplicate_rate_percent` are appended after `blastn_report`. To
diff against pre-2026-05-17 output, drop them first:

```
csvcut -C duplicate_pairs,duplicate_rate_percent,html_report \
    run_information_<batch>.csv | sort > new.csv
csvcut -C html_report old/run_information_<batch>.csv | sort > old.csv
diff old.csv new.csv   # should be empty
```

## Expected divergences

The following differences are intentional and should not be treated as
regressions:

- **fastp summary parsing path.** `virusHanter` scraped the fastp HTML
  with BeautifulSoup; `reportHanter` reads the fastp JSON directly. The
  values rendered in the "Read Summary from FASTP" tab are now formatted
  by `FastpProcessor` and may differ in unit prefixes and rounding from
  the strings the HTML scrape produced. The underlying read counts and
  Q20/Q30 fractions are the same.

- **Aggregation read-length column.** `virusHanter`'s aggregate parsed
  the read length out of the fastp HTML summary table; the refactored
  rule reads `summary.before_filtering.read1_mean_length` directly from
  the JSON. For paired-end reads with consistent length this matches; for
  ragged input the JSON value (mean) is the right answer.

- **Per-section error messages.** The refactored `ReportGenerator` wraps
  each section in its own try/except, so any per-section failure will
  produce a `Section '<name>' failed to build: ...` message rather than
  the generic `Failed to generate report` from the original.

## Sign-off

Document each parity run below as it is performed. Keep the refactored
pipeline in place as the default once at least two independent batches
have been validated.

| Date | Sample / batch | Operator | Result | Notes |
|------|----------------|----------|--------|-------|
|      |                |          |        |       |
