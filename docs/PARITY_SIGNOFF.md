# Parity sign-off procedure

How to validate that `virusHanter2 + reportHanter` reproduces the
original `virusHanter`'s `run_information_<batch>.csv`, and how to
record the result. Companion to
[`PARITY_NOTES.md`](PARITY_NOTES.md): that file catalogues the
intentional fixes and divergences; this file is the runnable checklist
whose result is pasted into the PARITY_NOTES sign-off table.

The sign-off table there is currently empty: parity is asserted by
design and guarded by unit tests, but has never been field-verified
against a real batch. This procedure closes that gap.

## What "parity" means here

`run_information_<batch>.csv` has 14 parity-locked columns. They are
not all compared the same way, because some are formatted strings, some
are numbers whose representation changed, and three are file blobs the
refactor deliberately altered. Group the columns by **how to compare
them**.

**Numeric (compare by value, not by formatting).** The underlying
number must match; the string representation may differ because the
refactor reads the fastp JSON instead of scraping the HTML:

- `number_reads`
- `read_len` (JSON `read1_mean_length`; matches for consistent
  paired-end input, may differ for ragged input -- a documented
  expected divergence)
- `mapped_to_human_percent`
- `kraken_virus_percent`
- `kaiju_virus_percent`

**Exact string (must match byte-for-byte).** Deterministic given the
same reference databases and the parity-recovery config below:

- `run_name`, `date`
- `top_virus_kaiju` (`"<taxon> (<reads>)"` joined by `||`, same Kaiju DB)
- `kaiju_report` (raw kaiju table blob; identical only when
  `TAXDUMP_NODES` is empty so no ICTV canonicalisation runs)

**Assembly-sensitive (review, do not auto-fail).** MEGAHIT can vary
slightly across versions/thread counts; a small benign difference is
not a parity failure -- investigate only large divergences:

- `number_of_contigs`, `top_contigs_blastn`

**Expected to differ (skip the comparison).** Different by design:

- `html_report` -- a different renderer (reportHanter); the blob never
  matches.
- `blastn_report` -- the per-assembler BLAST CSV carries an extra
  `assembler` column, so the blob differs structurally.

Trailing additive columns (`duplicate_pairs`, `assemblers_used`,
`host_removal_tool`, `<assembler>_n_contigs`, geNomad, provenance,
etc.) are dropped before diffing.

## Prerequisites

- **Linux.** CheckV 1.0.3 mis-reports "hmmsearch tasks failed" on macOS
  (see `rules/assembly.smk`); parity runs must be on Linux.
- **The original `virusHanter`** checked out and runnable, producing
  `run_information_<batch>.csv` for the batch.
- **The same reference databases for both runs.** This is the single
  most important prerequisite: `kraken_virus_percent`,
  `kaiju_virus_percent` and `top_virus_kaiju` are functions of the
  Kraken2 / Kaiju databases. Point both pipelines at the **identical**
  `KRAKEN_DB`, `KAIJU_DB`, `HUMAN_INDEX`, `BLASTN_DB` and
  `VIRUS_PARQUET`. If the databases differ, the data columns will
  differ for reasons unrelated to the refactor and the sign-off is
  meaningless.
- A representative batch: ideally >= 3 samples spanning a positive
  (clear viral hit), a host-heavy, and a near-negative sample.

## Parity-recovery config

Set these keys in the virusHanter2 config so every opt-in stage that
intentionally breaks parity is switched off:

```yaml
ASSEMBLERS: ["MEGAHIT"]            # single assembler (cols 9/10/14 from MEGAHIT only)
HOST_REMOVAL: "bwa"               # parity default; not hostile
DEDUPLICATE: "FALSE"             # no markdup filtering of reads
COVERAGE_SOURCES: ["KRAKEN"]      # pre-multi-source behaviour
COVERAGE_TOP_N: 20
COVERAGE_RANK_FILTER: []          # no rank filtering
COVERAGE_GENUS_WALKUP: "FALSE"
TAXDUMP_NODES: ""                # disables ICTV canonicalisation (keeps kaiju_report raw)
MULTIQC: "FALSE"
QUAST: "FALSE"
GENOMAD: "FALSE"
CLEAN: "FALSE"
# No SECONDARY_HOST_INDEX.
```

`CONTIG_LENGTH: 500` and `NUMBER_OF_PLOTS` should match whatever the
original run used.

## Run steps

1. Run the original `virusHanter` on the batch; keep its
   `run_information_<batch>.csv` as `original.csv`.
2. Run `virusHanter2` on the **same** batch with the parity-recovery
   config, against the **same** databases:

   ```bash
   conda activate virushanter
   cd virusHanter2
   snakemake --sdm conda --cores N --configfile config/config.parity.yaml
   ```

   Keep its `run_information_<batch>.csv` as `refactored.csv`.
3. Diff with the recipe below.

## Comparison recipe

Run this from any environment with pandas. It aligns on `sample_name`
and applies the comparison method appropriate to each column group.

```python
import re
import sys
import pandas as pd

INDEX = "sample_name"
NUMERIC_MUST = [
    "number_reads", "read_len", "mapped_to_human_percent",
    "kraken_virus_percent", "kaiju_virus_percent",
]
STRING_MUST = ["run_name", "date", "top_virus_kaiju", "kaiju_report"]
REVIEW = ["number_of_contigs", "top_contigs_blastn"]   # assembly-sensitive
SKIP = ["html_report", "blastn_report"]                # renderer / structure
TOL = 1e-6


def num(x):
    s = re.sub(r"[^0-9.eE+-]", "", str(x))
    try:
        return float(s)
    except ValueError:
        return float("nan")


orig = pd.read_csv(sys.argv[1]).set_index(INDEX).sort_index()
new = pd.read_csv(sys.argv[2]).set_index(INDEX).sort_index()

if list(orig.index) != list(new.index):
    print("FAIL: sample sets differ")
    print("  original:  ", list(orig.index))
    print("  refactored:", list(new.index))
    sys.exit(1)

fails, reviews = [], []
for col in NUMERIC_MUST:
    for s in orig.index:
        o, n = num(orig.at[s, col]), num(new.at[s, col])
        if not abs(o - n) <= TOL * max(1.0, abs(o)):
            fails.append((s, col, orig.at[s, col], new.at[s, col]))
for col in STRING_MUST:
    for s in orig.index:
        if str(orig.at[s, col]) != str(new.at[s, col]):
            fails.append((s, col, orig.at[s, col], new.at[s, col]))
for col in REVIEW:
    for s in orig.index:
        if str(orig.at[s, col]) != str(new.at[s, col]):
            reviews.append((s, col, orig.at[s, col], new.at[s, col]))

print(f"FAILS (must match): {len(fails)}")
for s, col, o, n in fails:
    print(f"  FAIL {s} | {col}: {o!r} -> {n!r}")
print(f"REVIEW (assembly-sensitive): {len(reviews)}")
for s, col, o, n in reviews:
    print(f"  REVIEW {s} | {col}: {o!r} -> {n!r}")
print(f"SKIPPED (different by design): {', '.join(SKIP)}")

sys.exit(0 if not fails else 2)
```

```bash
python parity_diff.py original.csv refactored.csv
```

**Pass condition:** **0 FAILS** (exit 0). Any REVIEW rows are small and
explained (MEGAHIT non-determinism / version). SKIPPED columns are not
compared.

## Checklist

- [ ] Run on Linux.
- [ ] Both pipelines pointed at the identical Kraken2, Kaiju, HUMAN,
      BLAST and VIRUS_PARQUET databases (record the paths / build dates).
- [ ] Parity-recovery config applied (`ASSEMBLERS: ["MEGAHIT"]`,
      `TAXDUMP_NODES: ""`, multi-source / rank-filter / dedup off, etc.).
- [ ] Original `virusHanter` run completed -> `original.csv`.
- [ ] `virusHanter2` run completed -> `refactored.csv`.
- [ ] Same sample set in both CSVs.
- [ ] `parity_diff.py` reports **0 FAILS**.
- [ ] Any REVIEW (assembly-sensitive) differences reviewed and explained.
- [ ] Result recorded in the PARITY_NOTES sign-off table.
- [ ] Repeat on a **second independent batch** before declaring the
      refactored pipeline the validated default (PARITY_NOTES requires
      two batches).

## Recording the result

Append a row to the sign-off table in
[`PARITY_NOTES.md`](PARITY_NOTES.md) (the empty `## Sign-off` table at
the end), for example:

```
| 2026-06-26 | 251015_..._DRRKK (3 samples) | <operator> | PASS (0 FAILS; contigs +/-1 on sample 138) | DBs: k2_viral_20260517, kaiju_refseq_viral_20260517 |
```

State the database snapshot used and any REVIEW note. Two PASS rows on
independent batches close the production parity gate.
