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
`report_html_blob` column (hex-encoded HTML, which carries the same
timestamp/build-id non-determinism described above). Diff the two CSVs
with that column dropped, e.g.:

```
csvcut -C report_html_blob run_information_<batch>.csv | sort > a.csv
csvcut -C report_html_blob run_information_<batch>.csv | sort > b.csv
diff a.csv b.csv
```

Numeric columns (read counts, percentages, `number_of_contigs`) should
match exactly. String columns (`top_contigs_blastn`, `top_virus_kaiju`)
should match exactly unless one of the documented divergences below
applies.

## Intentional fixes from the original

The following number is *not* a parity divergence — it is a fix to a
pre-existing aggregation bug in the original `virusHanter`:

- **`kraken_virus_percent`**. The original `virusHanter` ran
  `kraken_df.loc[kraken_df['domain'] == 'Viruses', 'percent'].sum()`,
  which sums every row in the Kraken report whose `domain` column is
  "Viruses" — i.e. the Domain (D) row plus every species (S) row
  underneath it. Since the Domain row's percent already accounts for
  every clade beneath it, this double-counts the viral fraction (often
  by 2-3x, depending on how deep the species hierarchy reaches).
  `virusHanter2` instead reads the single Domain-level row directly.
  The new value is the correct viral percentage; an old `run_information_*.csv`
  with this column will show a strictly larger (often 2-3x) value.

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

- **Aggregation column order.** Column names are preserved; column order
  is determined by Python `dict` insertion order, which is stable but may
  differ from the original.

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
