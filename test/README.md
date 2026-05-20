# virusHanter2 smoke testing

Three test tiers live here.

## Fixture overview

The smoke pipeline runs against three synthetic samples, each
carrying one of three pseudo-random 5 kb "viruses":

| Sample    | Carries virus | Taxid  | Virus reads / coverage |
|-----------|---------------|--------|------------------------|
| sample1_R | alpha         | 100001 | 800 pairs (~50x)       |
| sample2_R | beta          | 100002 | 600 pairs (~36x)       |
| sample3_R | gamma         | 100003 | 400 pairs (~24x)       |

The three viruses are non-repetitive and seeded independently so MEGAHIT
assembles each one to a single contig and the classifiers (Kraken2,
Kaiju, BLASTN) can tell them apart. The fixture exercises both the
multi-sample paths (per-sample reports, run-level aggregation,
per-virus aggregation) and the substring-attribution path that links
contigs to taxids in `scripts/per_virus_metrics.py`.

## Dry-run (no databases needed)

```
./test/run_smoke.sh
```

Runs `snakemake --lint` and `snakemake -n --sdm conda` against
`test/config.test.yaml`. No tools are invoked, so nothing under
`test/mini_db/` has to exist. This is the cheapest regression check
for DAG construction and the rule imports.

The fixture FASTQs (`sample1_R{1,2}.fastq.gz`, `sample2_R{1,2}.fastq.gz`,
`sample3_R{1,2}.fastq.gz`) are not committed; a fresh checkout has
none of them on disk. Run `./test/build_fixtures.sh` (or
`./test/run_smoke.sh --build`) to synthesise them alongside the mock
databases.

## Partial smoke (synthetic mini-DBs, no CheckV)

```
./test/run_smoke.sh --full     # auto-degrades when CheckV is stubbed
```

`build_fixtures.sh` synthesises everything in `test/mini_db/` except a
real CheckV database. The mock databases all carry the three viruses
so the Kraken2 and Kaiju mini-DBs cover taxids 100001-100003. Tools
required on `$PATH` to build the fixtures:

| Tool | Used for |
|---|---|
| `python` + `pandas` + `pyarrow` | FASTQ synthesis, `virus.parquet` |
| `bwa` | host BWA index |
| `kraken2-build` | Kraken2 mini-DB (three taxids) |
| `kaiju-mkbwt` + `kaiju-mkfmi` | Kaiju mini-DB (three proteins) |
| `makeblastdb` | BLAST nt mini-DB (three records) |

When `test/mini_db/checkv` is missing or contains only a `.stub`
sentinel, the smoke runs through `--until blastn mosdepth_kraken_hits
kaiju_to_table` for all three samples and asserts the BLASTN, KAIJU,
KRAKEN and MOSDEPTH outputs exist for each. The HTML reports and the
run-information aggregation are skipped (they depend on
`merge_checkv_blastn`).

## Full smoke (real CheckV database)

Provide a real CheckV database at `test/mini_db/checkv/` (e.g. by
downloading it with `checkv download_database test/mini_db/`) or
symlink one in:

```
mkdir -p test/mini_db/checkv
ln -s /path/to/checkv-db-v1.5/genome_db test/mini_db/checkv/genome_db
ln -s /path/to/checkv-db-v1.5/hmm_db    test/mini_db/checkv/hmm_db
```

`./test/run_smoke.sh --full` then runs the complete pipeline,
including `generate_report` and the run-level aggregation, and
asserts the per-sample HTML and the run-information CSV exist for
all three samples.

### macOS AppleDouble gotcha

If the CheckV DB sits on an external Mac-mounted volume (HFS+ or APFS
formatted by Finder), macOS may sprinkle AppleDouble metadata files
(`._<name>.hmm`) into `hmm_db/checkv_hmms/`. CheckV's hmmsearch
driver picks these up as inputs and reports "80 hmmsearch tasks
failed" even when every real HMM ran cleanly. Strip them once with:

```
find /path/to/checkv-db-v1.5 -name '._*' -delete
```

See `docs/REFERENCE_DBS.md` for a longer note. The diagnosis was
worked out in 2026-05-20 after a previous (mistaken) attribution to
a `mp.Pool` / `hmmsearch` interaction.
