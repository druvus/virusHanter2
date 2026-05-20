# virusHanter2 smoke testing

Three test tiers live here.

## Dry-run (no databases needed)

```
./test/run_smoke.sh
```

Runs `snakemake --lint` and `snakemake -n --sdm conda` against
`test/config.test.yaml`. No tools are invoked, so nothing under
`test/mini_db/` has to exist. This is what CI runs on every commit and is
the cheapest regression check for DAG construction and the rule imports.

The fixture FASTQs (`sample1_R{1,2}.fastq.gz`, `sample2_R{1,2}.fastq.gz`,
`sample3_R{1,2}.fastq.gz`) are not committed; a fresh checkout has none
of them on disk. Run `./test/build_fixtures.sh` (or
`./test/run_smoke.sh --build`) to synthesise them alongside the mock
databases. Each sample uses a distinct random seed and read mix so the
smoke exercises multi-sample paths (per-sample reports plus run-level
aggregation).

## Partial smoke (synthetic mini-DBs, no CheckV)

```
./test/run_smoke.sh --full     # auto-degrades when CheckV is stubbed
```

`build_fixtures.sh` synthesises everything in `test/mini_db/` except a real
CheckV database. Tools required on `$PATH` to build the fixtures:

| Tool | Used for |
|---|---|
| `python` + `pandas` + `pyarrow` | FASTQ synthesis, `virus.parquet` |
| `bwa` | host BWA index |
| `kraken2-build` | Kraken2 mini-DB |
| `kaiju-mkbwt` + `kaiju-mkfmi` | Kaiju mini-DB |
| `makeblastdb` | BLAST nt mini-DB |

When `test/mini_db/checkv` contains the `.stub` sentinel, the smoke runs
`snakemake --until blastn` and asserts the BLASTN output exists. The HTML
report and run-info aggregation are skipped (they depend on
`merge_checkv_blastn`).

## Full smoke (real CheckV database)

Provide a real CheckV database at `test/mini_db/checkv/` (e.g. by
downloading it once with `checkv download_database test/mini_db/`) and
remove the `.stub` sentinel. `./test/run_smoke.sh --full` will then run
the complete pipeline including `generate_report` and assert the per-sample
HTML and the run-information CSV exist.
