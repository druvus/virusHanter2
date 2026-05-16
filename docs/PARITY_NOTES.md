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
