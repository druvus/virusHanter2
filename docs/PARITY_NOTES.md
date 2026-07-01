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
   snakemake --sdm conda --cores <N> --configfile config/config.yaml
   ```

3. Compare outputs for each sample.

### Per-sample HTML

The HTML reports are not byte-comparable (Panel embeds build-time
timestamps and randomised DOM ids). Compare instead by opening both files
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
  is 500. (The `PLOT_THRESHOLD` knob from the original `virusHanter`
  was retired together with the `bam2plot` rule; see the 2026-05-19
  section below.)
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
  emits. Since the 2026-05-19 retirement of `bam2plot`, the
  `regions.bed.gz` is also what drives the interactive coverage
  traces in the per-sample HTML report.
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

## Additive 2026-05-19 viralrecon-gap additions

Two further opt-in rules were added after a feature gap analysis
against nf-core/viralrecon. Both are off by default so the parity
invariant continues to hold for existing runs; flip the flag to opt
in.

- **`DEDUPLICATE: "TRUE"` flag**
  (`rules/pre_processing.smk`). When set, `markdup_human` writes a
  marked BAM and `remove_host` reads it with `samtools view -F 1024`
  so PCR duplicates do not propagate into the host-removed FASTQs
  that feed MEGAHIT and the BWA-to-Kraken-hits coverage step. With
  the flag off (default), `remove_host` continues to consume the
  un-marked `bwa_human` BAM and outputs are byte-identical to before.
  When enabling on a production sample, re-run an existing sample
  with the flag both off and on and compare the resulting Kraken /
  Kaiju top-hit percentages and mosdepth mean coverages — the
  magnitude of the shift is library-dependent and should be recorded
  here for the panels that route through this option.
- **`QUAST: "TRUE"` flag** (`rules/assembly.smk`). When set,
  `quast_megahit` runs against each sample's raw MEGAHIT contigs and
  writes `QUAST/report.tsv` per sample. The report directory is
  picked up automatically by MultiQC and also depended on as an
  explicit input to keep the dependency graph honest. With the flag
  off (default), no QUAST process runs and there is no change to the
  rest of the workflow.

Both additions are pure additions to the workflow graph: the
parity-locked columns of `run_information_<batch>.csv` are not
touched in either case, and the per-sample HTML report is unchanged
unless and until reportHanter is taught to surface QUAST metrics
(separate follow-up).

## 2026-05-20: Apple Silicon MEGAHIT k=21 SIGSEGV — `--k-min 27`

The smoke pipeline's apparent "CheckV failure on Apple Silicon" was
re-diagnosed. `megahit_core_no_hw_accel count -k 21` SIGSEGVs on
small inputs on osx-arm64 (the same kernel-of-bug as the `_popcnt`
variant). The MEGAHIT rule's safety net then wrote `DUMMY_CONTIG`,
which Pilon polished into a 200 bp repeat. Prodigal-gv (called by
CheckV) found no ORFs in the dummy and emitted an empty
`proteins.faa`; CheckV's hmmsearch then errored on the empty input
with "Sequence file is empty or misformatted" for all 80 tasks,
which the earlier session mistook for an `mp.Pool` internal failure
(observations #1987, #1990, #1992).

The fix in `rules/assembly.smk` adds `--k-min 27` on Apple Silicon,
alongside the existing `--no-hw-accel` flag and the 2-thread cap.
That avoids the buggy k=21 path entirely. Verified: MEGAHIT runs
the smoke input to a real 4906 bp contig in ~0.2 s; CheckV then
produces a real `contamination.tsv` in ~7 s.

Production Linux runs are unaffected — both `_popcnt` and
`_no_hw_accel` binaries handle k=21 correctly there. The `--k-min 27`
flag is gated on the platform check, so it imposes no
assembly-quality concession on Linux.

## 2026-05-21: Rank filter + genus walk-up on the coverage reference set

`bwa_align_to_kraken_hits` now consults an NCBI taxdump
(`nodes.dmp`) when `TAXDUMP_NODES` is configured. Two new
behaviours are gated on the taxdump:

1. **Rank filter** (`COVERAGE_RANK_FILTER`, default
   `[realm, kingdom, subkingdom, phylum, subphylum, class,
   subclass, order, suborder, family, subfamily]`). Classifier
   hits at these higher ranks are dropped silently before they
   enter the coverage union. Without the filter, NCBI taxonomy
   propagation rows like `Viruses` (kingdom), `Cardeaviricetes`
   (class), `Anelloviridae` (family) and `Herpesvirales` (order)
   flood `unmapped_taxids.tsv` because no per-taxid sequence
   exists for ranks above genus.

2. **Genus walk-up** (`COVERAGE_GENUS_WALKUP`, default `TRUE`).
   When a species/strain taxid is absent from `VIRUS_PARQUET`
   the rule walks its parent chain via `nodes.dmp` to the first
   ancestor at rank `genus`. If a parquet reference exists for
   that genus the rule substitutes it and tags the source in
   `virus_names` with `->genus` (e.g. `kraken->genus`). Two
   strain-level Kraken hits that share a genus collapse to a
   single genus reference, which is the desired behaviour: the
   coverage panel surfaces "any Alphatorquevirus" coverage even
   when the specific strains the classifier called are not in
   the parquet.

This is a **deliberate parity break**:

- The `virus_names` sidecar's `sources` column can now carry
  `->genus`-suffixed tags. Older sidecars that pre-date this
  change still render because the suffix is only an additional
  string value.
- The mosdepth chrom set may include genus-level references
  that were never in classifier output directly. Coverage tab
  labels like `NC_038338.1 — Alphatorquevirus [kraken->genus]`
  signal the substitution.

To recover today's behaviour, set `COVERAGE_RANK_FILTER: []`
and `COVERAGE_GENUS_WALKUP: "FALSE"`, or leave `TAXDUMP_NODES`
empty (the rule then logs a single warning and degrades to the
multi-source union without rank filter or walk-up).

The build-time parquet rebuild now writes two extra columns
(`rank`, `genus_taxid`) when the refresh workflow points at a
`nodes.dmp`. The columns are purely additive — every existing
consumer (`parquet_accession_to_taxid`,
`parquet_refs_by_taxid`, `pick_longest_per_taxid`) accesses the
DataFrame by column name and ignores extras.

## 2026-05-21: Multi-source coverage reference set

`bwa_align_to_kraken_hits` previously picked the Kraken2 top-20
viral taxa, intersected them with `VIRUS_PARQUET`, and aligned the
host-removed reads to the resulting FASTA. The mosdepth chrom set
therefore reflected only Kraken's evidence, and only that subset
of Kraken's taxa for which `VIRUS_PARQUET` had a reference. Kaiju
hits never entered the coverage step. BLASTN hits from the
assembled contigs never entered the coverage step either, so a
sample where the assemblers recovered a herpesvirus contig still
got no coverage trace for that virus when Kraken's herpesvirus
percent fell outside the hard-coded top-20.

The rule's `run:` block has been rewritten to take the union of
the configured `COVERAGE_SOURCES` (default
`["KRAKEN", "KAIJU", "BLAST"]`):

- KRAKEN: as before, top-`COVERAGE_TOP_N` by percent (default 20).
- KAIJU: top-N by percent from `kaiju_to_table`, filtered against
  `VIRUS_PARQUET`'s tax_id set so non-viral RefSeq hits do not
  enter.
- BLAST: every per-assembler merged CSV; accessions are resolved
  to taxids via the new `parquet_accession_to_taxid` helper in
  `scripts/functions.py`.

The output sidecar `virus_names` gains a trailing `sources`
column (`kraken`, `kaiju`, `blast`, or semicolon-delimited
combinations). A new `unmapped_taxids.tsv` sidecar per sample
lists classified taxids that had no reference in
`VIRUS_PARQUET`, so the reviewer can see which evidence was
dropped and which families (e.g. specific Anelloviridae strains)
need a parquet rebuild.

This is a **deliberate parity break**:

- The mosdepth chrom set is no longer Kraken-top-20 only; it is
  the union over `COVERAGE_SOURCES`. The Alignment Coverage panel
  in the per-sample report shows more tabs (typically 20–50
  rather than 5–15).
- The `virus_names` sidecar grows from three columns to four;
  reportHanter v0.5.1 reads the optional `sources` column
  back-compat with three-column sidecars from older runs.

To recover byte-identical behaviour with the pre-this-change
runs, set `COVERAGE_SOURCES: ["KRAKEN"]` and
`COVERAGE_TOP_N: 20`. The `virus_names` file then still carries
the new column but every row's `sources` value is `kraken`, so
the column can be dropped before diffing.

## 2026-05-21: Multi-assembler mode (MEGAHIT + metaSPAdes + rnaviralSPAdes)

`ASSEMBLERS: ["MEGAHIT", "metaSPAdes", "rnaviralSPAdes"]` is now the
default (all three; the deprecated `SPAdes` alias is rejected at
workflow load). Every contig-producing rule (Pilon, BLASTN, CheckV,
geNomad and QUAST when enabled) runs once per (sample, assembler) and
lands under
`{sample}/{assembler}/...`. The report's "Classification of Contigs"
tab carries an `assembler` column and the BLAST headline bar chart
splits per assembler. The per-virus CSV gains trailing
`{assembler}_contigs` columns; the run-info CSV gains
`assemblers_used` plus `{assembler}_n_contigs` / `{assembler}_n50`.

This is a **deliberate parity break**:

- `number_of_contigs` in `run_information_<batch>.csv` now sums
  contigs across all active assemblers, so the value will differ
  from a MEGAHIT-only baseline.
- `top_contigs_blastn` selects from the union; the top five may
  include rows from any of the active assemblers.
- The HTML report's contig table now has an extra leading column.

To recover byte-identical parity with the original `virusHanter`,
set `ASSEMBLERS: ["MEGAHIT"]` in the config. Every other rule then
behaves exactly as before: paths still carry the `{assembler}`
wildcard but the wildcard takes a single value, so the BLAST /
CheckV merged CSVs sit at `{sample}/MEGAHIT/BLASTN/...` rather than
the original `{sample}/BLASTN/...`. The CSV columns added in this
session are all trailing and can be dropped before diffing.

Apple Silicon: bioconda has `osx-arm64` builds for `spades>=3.15.5`
so metaSPAdes runs alongside MEGAHIT on a Mac without
`CONDA_SUBDIR=osx-64`. The rule mirrors MEGAHIT's "dummy contig on
failure" fallback so a SPAdes refusal (it imposes a per-library
minimum) does not tear the DAG down.

## Additive: provenance (databases + application versions)

The run records which reference databases and application versions
produced it (see [PROVENANCE.md](PROVENANCE.md)). This is **purely
additive**:

- `run_information_<batch>.csv` gains two trailing columns after the
  existing provenance block: `databases_build_identity` (one
  `KEY=<identity>` per DB, preferring a robust build stamp over an
  mtime) and `tool_versions` (the conda-resolved headline tool
  versions, e.g. `fastp=0.24.0;kraken2=2.1.3`). Both are blank on a
  legacy run, so a column-dropped diff stays clean. The parity-locked
  first 14 columns are untouched.
- Two new run-level outputs are produced, `software_versions.tsv` and
  `run_provenance_<batch>.json`; neither feeds another rule nor alters
  an existing one.
- The HTML report gains a `Provenance` tab. The rendered HTML is
  already an "expected to differ" column in parity, so this does not
  affect the locked schema.

Drop the new trailing columns before diffing:

```
csvcut -C databases_build_identity,tool_versions,html_report \
    run_information_<batch>.csv | sort > new.csv
```

## 2026-05-19: `bam2plot` retired

The `bam2plot` rule that produced `COVERAGE_PLOTS/*.svg` and the
matching `PLOT_THRESHOLD` config key have been removed. Since
`mosdepth_kraken_hits` now feeds an interactive coverage trace per
reference into the report directly (via the per-sample
`regions.bed.gz`), the SVGs were never read and added a heavy
dependency that segfaulted on Apple Silicon. Practical impact:

- `envs/bam2plot.yaml` deleted; the smoke runner no longer carries
  Apple-Silicon-specific workarounds.
- `--coverage_folder` is no longer a flag of the `reporthanter` CLI
  or a parameter of `create_report` / `ReportGenerator`. The CLI now
  requires `--mosdepth_regions`.
- `COVERAGE_PLOTS/` directories from previous runs can be removed
  freely; nothing in the current workflow reads them.
- Parity-locked columns of `run_information_<batch>.csv` are
  unaffected because none of them sourced data from `bam2plot`.

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

Follow [`PARITY_SIGNOFF.md`](PARITY_SIGNOFF.md) for the runnable
procedure (parity-recovery config, run steps, the `parity_diff.py`
comparison recipe and the pass condition). Document each parity run
below as it is performed. Keep the refactored pipeline in place as the
default once at least two independent batches have been validated.

| Date | Sample / batch | Operator | Result | Notes |
|------|----------------|----------|--------|-------|
|      |                |          |        |       |
