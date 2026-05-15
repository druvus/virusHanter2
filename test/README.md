# virusHanter2 smoke testing

Two test tiers live here.

## Dry-run (no databases needed)

```
./test/run_smoke.sh
```

Runs `snakemake --lint` and `snakemake -n --use-conda` against
`test/config.test.yaml` using the committed `test_R[12].fastq.gz` reads. The
dry-run does not invoke any tool, so the mock databases under
`test/mini_db/` do not have to exist on disk for it to pass. This is what CI
runs on every commit and is the cheapest regression check for the DAG and
the rule imports.

### Caveat about the committed fixture FASTQs

`test_R[12].fastq.gz` ship as empty placeholders. The sample-discovery helper
`paired_reads()` only matches files whose extension is `.fq`, `.fastq`, `.fa`,
`.fasta`, or `.fna` (a behavior inherited from `virusHanter`), so as long as
the placeholders are gzipped the dry-run reports an empty per-sample plan and
only schedules the aggregate rule. To exercise the per-sample chain
end-to-end, replace the placeholders with real un-gzipped reads (or add a
small synthetic pair).

## Full smoke (mini reference databases required)

```
./test/run_smoke.sh --full
```

Requires the following minimal databases under `test/mini_db/`:

- `human/` — BWA index of a single short reference chromosome; prefix is `human`.
- `kaiju/` — directory containing one `.fmi` file plus `names.dmp` and `nodes.dmp`.
- `kraken/` — Kraken2 mini-DB (built with `kraken2-build --special viral` or similar).
- `blast/` — BLAST nucleotide database with prefix `viral` (use `makeblastdb -in viruses.fasta -dbtype nucl -out viral`).
- `checkv/` — CheckV mini-DB. Construct manually or symlink to a real one.
- `virus.parquet` — Parquet of `name`, `sequence`, `tax_id` columns for a handful of viruses.

None of these files are tracked in git; build or copy them in locally before
running `--full`. The full run produces `test/results/test/test/REPORT/test.html`
and `test/results/test/run_information_test.csv` and the script asserts both
exist.
