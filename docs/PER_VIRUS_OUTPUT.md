# Per-(sample, virus) CSV and multi-run master files

The collaborator-facing output is `per_virus_<batch>.csv` written at
the top of each batch's result folder. One row per detected Kraken
viral taxid, per sample, capped at `NUMBER_OF_PLOTS` (config;
default 10).

`bwa_align_to_kraken_hits` itself builds a BWA index over a wider
Kraken top-20 (hardcoded, to match the original `virusHanter`); the
narrower `NUMBER_OF_PLOTS` cap is applied downstream by
`per_virus_metrics`. Coverage traces for the same references are
rendered in the per-sample HTML report from the mosdepth
`regions.bed.gz` produced by `mosdepth_kraken_hits`.

## Per-sample, per-virus schema

| Column | Source / formula |
|---|---|
| `run_name` | Batch folder name (the basename of `SAMPLES`). |
| `sample_name` | Per-sample folder name. |
| `date` | `run_name.split("_")[0]`. |
| `virus_name_kraken` | Kraken2 report `name` for this taxon. |
| `virus_taxid` | Kraken2 `taxonomy_id`. Row key together with `sample_name`. |
| `virus_name_kaiju` | Kaiju TSV `taxon_name` joined by `taxon_id`. Empty if no match. |
| `contigs` | Count of BLASTN-merged contigs attributed to this taxid (accession first, name substring fallback). |
| `virus_reads_kraken2` | Kraken2 `count_clades` for this row. |
| `other_reads` | `non_human_reads - all_viral_kraken_reads`. |
| `total_reads` | fastp `summary.before_filtering.total_reads`. |
| `human_reads` | flagstat `with itself and mate mapped`. |
| `human_reads_percent` | `100 * human_reads / total_reads`. |
| `non_human_reads` | `total_reads - human_reads`. |
| `non_human_reads_percent` | `100 - human_reads_percent`. |
| `note` | Empty for clean runs; carries `"MEGAHIT assembly failed; dummy contig only"` when every BLASTN contig for the sample is the `DUMMY_CONTIG` fallback so silent assembly failures are visible. |
| `specific_virus_rpm` | `virus_reads_kraken2 / total_reads * 1e6`. |
| `all_virus_rpm` | `all_viral_kraken_reads / total_reads * 1e6`. `all_viral_kraken_reads` is the Kraken Domain "Viruses" row's `count_clades` (already accounts for descendant clades). |
| `Completeness (% >5X)` | Percent of reference bases with mosdepth depth >= 5 across all references for this taxid. 0-100. |
| `bases_above_5x` | Raw base count >= 5x across all references for this taxid. |
| `mean_coverage` | Weighted mean depth: `sum(bases_aligned) / sum(reference_length)` across references. |

### Trailing additive columns

Appended after the parity-locked columns above; tolerated by all
existing consumers because they read the CSV by column name.

| Column | Source / formula |
|---|---|
| `<assembler>_contigs` | One column per active entry in `config[ASSEMBLERS]`. Count of BLASTN-merged contigs attributed to this taxid **for that assembler**. Sums to `contigs` across assemblers. |
| `genomad_viral_contigs` | When `GENOMAD: "TRUE"`: count of attributed contigs geNomad called viral. |
| `genomad_max_virus_score` | When `GENOMAD: "TRUE"`: highest geNomad `virus_score` among the attributed contigs. |

A virus that is in the Kraken top-N but has no aligned reads ends up
with `bases_above_5x = 0`, `Completeness (% >5X) = 0`, `mean_coverage = 0`.
The "Viruses" Domain row itself is one of the top-N candidates and
has no reference in the BWA index; its row carries the Kraken stats
but blank coverage.

## Multi-reference taxids

`VIRUS_PARQUET` can have several reference sequences per taxid. The
per-virus row aggregates across them:

- `bases_above_5x` = sum across references.
- `mean_coverage` = `sum(bases_aligned) / sum(reference_length)`.
- `Completeness (% >5X)` = `100 * bases_above_5x / sum(reference_length)`.

## Per-batch concatenation

`scripts/aggregate_per_virus.py` runs as the `aggregate_per_virus`
Snakemake rule and produces `per_virus_<batch>.csv` in the result
folder. Schema and column order match the per-sample file. Run on
demand via:

```
snakemake --sdm conda --cores N --configfile config/config.production.yaml \
    --until aggregate_per_virus
```

(`per_virus_<batch>.csv` is also a `rule all` target, so a full
workflow run produces it without a `--until`.)

## Multi-run master files

`scripts/merge_runs.py` is a standalone CLI. It does not call
Snakemake; it only reads per-batch CSVs and writes per-multi-run
master CSVs.

```
python scripts/merge_runs.py \
    --result-folder /path/to/RESULTS/<batch1> \
    --result-folder /path/to/RESULTS/<batch2> \
    --result-folder /path/to/RESULTS/<batch3> \
    --out-dir /path/to/master/
# writes:
#   /path/to/master/master_per_sample.csv
#   /path/to/master/master_per_virus.csv
```

`master_per_sample.csv` concatenates every batch's
`run_information_<batch>.csv`. `master_per_virus.csv` concatenates
every batch's `per_virus_<batch>.csv`. Each input row carries its own
`run_name`, so the master files are self-describing.

## Where the numbers come from

- Kraken2 report -> `wrangle_kraken` CSV -> per-virus virus name,
  taxid, `virus_reads_kraken2`, and (Domain row only) the
  `all_viral_kraken_reads` numerator.
- Kaiju TSV -> `virus_name_kaiju` joined by `taxon_id`.
- BLASTN merged CSV (with CheckV inner join) -> `contigs` via either
  parquet-accession match or first-token substring against the
  Kraken taxon name. With `config[ASSEMBLERS]` carrying more than
  one entry the script reads every per-assembler merged CSV and
  partitions the matched contigs across the trailing
  `<assembler>_contigs` columns; the parity-existing `contigs`
  column sums across assemblers.
- mosdepth `summary.txt` -> per-reference length, total bases, mean.
- mosdepth `thresholds.bed.gz` (from `--thresholds 1,5,10`) ->
  per-reference `bases_above_5x` (sum of the 5X column per chrom).
- fastp JSON -> `total_reads`.
- Host flagstat -> `human_reads` (paired with mate mapped count).

The per-sample summary `run_information_<batch>.csv` is unchanged
and untouched by this addition.
