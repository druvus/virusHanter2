# Tutorial: rebuild the classification databases from one snapshot

The classification stack in `virusHanter2` is held together by the
NCBI taxonomy ID:

- **`VIRUS_PARQUET`** carries one nucleotide reference per viral
  tax_id and drives the BWA + mosdepth coverage step.
- **`KAIJU_DB`** is a Burrows–Wheeler protein index whose
  sequences are tagged with NCBI tax_ids.
- **`KRAKEN_DB`** is a *k*-mer-based DNA index whose leaves are
  NCBI tax_ids.
- **`BLASTN_DB`** is a nucleotide BLAST database (alias) whose
  hits are mapped back to tax_ids via the `taxdb` files NCBI
  publishes alongside.
- **`TAXDUMP_NODES`** (optional) is the NCBI `nodes.dmp` the
  coverage rule uses for the rank filter and the genus walk-up.

When these databases are built from different NCBI snapshots, a
read can be classified to a taxon by Kraken2 that has no
representative reference in the parquet, or Kaiju can hit a
protein under a taxon that Kraken2 does not know exists yet. The
existing rank filter + genus walk-up shipped in the pipeline
absorbs the residual asymmetry at runtime, but starting from a
single coordinated snapshot keeps the diagnostic load light.

The refresh workflow at
[`refresh/refresh_virus_parquet.smk`](../refresh/refresh_virus_parquet.smk)
rebuilds the parquet, the Kaiju FM-index and the Kraken2 viral
DB from the same NCBI viral RefSeq snapshot, downloads the
matching taxdump (`nodes.dmp` + `names.dmp`), and emits an
overlap-with-Kraken2 sidecar so the operator can see at a glance
which tax_ids are covered by which classifier. The workflow also
drives `update_blastdb.pl` so the viral BLAST DB tarballs share
the same snapshot as the parquet, the Kaiju FMI and the Kraken2
DB; no separate manual step is required. Building Kraken2 from
the same RefSeq pull closes the recurring gap where the publicly
hosted `k2_viral_*` snapshots occasionally omit individual
genomes (e.g. the Feb 2026 snapshot missed HSV-2 / NC_001798
while the matching parquet and Kaiju FMI both included it).

## What you need before you start

- A LaCie-class external drive (or any volume with ~40 GB free
  for downloads and ~5 GB for the published outputs).
- The `virushanter` conda env active (or any env with
  `snakemake-minimal=9.23.*` plus `pandas`, `pyarrow`, `pyfastx`,
  `requests`-stack — the refresh workflow materialises its own
  per-rule conda env from `envs/refresh.yaml` regardless).
- ~1 hour of wall time on a residential connection: the longest
  pole is the ~11 GB `prot.accession2taxid.gz` download.

## One-shot refresh

```bash
conda activate virushanter
cd virusHanter2

# First time only: copy the example config and edit paths.
cp refresh/config.yaml refresh/config.local.yaml
$EDITOR refresh/config.local.yaml

# Run. Re-running with the same config is idempotent thanks to
# Snakemake's mtime tracking; pass --rerun-incomplete if a previous
# attempt was interrupted.
snakemake -s refresh/refresh_virus_parquet.smk \
    --configfile refresh/config.local.yaml --cores 4 \
    --sdm conda
```

`refresh/config.local.yaml` carries four key paths:

```yaml
OUTPUT_PARQUET: "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet"
DOWNLOAD_DIR:   "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/_refresh_workdir"
KRAKEN_DB_FOR_COMPARE: "/Volumes/LaCie/REGIONEN/ref_dbs/KRAKEN_DB/k2_viral_20260226"
# Source URLs default to the canonical NCBI FTP paths; override
# only when you want a pinned snapshot.
```

## What the workflow does

The DAG follows the `download_*` / `extract_*` / `build_*` /
`publish_*` rule pattern. From a clean working directory the
sequence is:

1. **Download** the multi-part viral RefSeq genomic FASTA
   (`viral.*.1.genomic.fna.gz`), the matching protein FASTA
   (`viral.*.protein.faa.gz`), the nucleotide and protein
   accession2taxid mappings, and the NCBI taxdump tarball. The
   rules use `curl` (not `wget`) so the `utime()` LaCie filesystem
   quirk does not abort the chain.
2. **Decompress** the parts into single concatenated FASTAs and
   extract `nodes.dmp` + `names.dmp` from `taxdump.tar.gz`.
3. **Build the parquet** via
   [`scripts/build_virus_parquet.py`](../scripts/build_virus_parquet.py)
   with `--source refseq --one-rep-per-taxid --taxdump-nodes nodes.dmp`.
   Each row carries `name`, `sequence`, `tax_id`, `rank`,
   `genus_taxid`.
4. **Build the Kaiju FMI index** from the protein FASTA. Headers
   are rewritten to the bare-tax_id format Kaiju's `kaiju-mkbwt`
   actually expects (see the dedicated header rewriter at
   [`scripts/reformat_kaiju_headers.py`](../scripts/reformat_kaiju_headers.py)),
   then `kaiju-mkbwt` and `kaiju-mkfmi` build the index. The
   resulting `kaiju_refseq_viral.fmi` is published next to the
   parquet alongside `nodes.dmp` and `names.dmp`.
5. **Build the Kraken2 viral DB** from the same nucleotide
   FASTA. The rule seeds Kraken2's `taxonomy/` folder with the
   already-downloaded `nodes.dmp`, `names.dmp` and decompressed
   `nucl_gb.accession2taxid`, then runs `kraken2-build
   --add-to-library` followed by `kraken2-build --build` to hash
   the k-mers. The built DB (`hash.k2d`, `taxo.k2d`, `opts.k2d`,
   `seqid2taxid.map`, `inspect.txt`) is published to
   `kraken2_refseq_viral/` next to the parquet so the main
   pipeline's `KRAKEN_DB` key can point at a snapshot that
   matches the Kaiju FMI's taxid universe.
6. **Publish** the `nodes.dmp` to a stable path next to the
   parquet so the main pipeline's `TAXDUMP_NODES` config key can
   point at it without re-extracting the tar.
7. **Compare with Kraken2** via
   [`scripts/compare_parquet_kraken2.py`](../scripts/compare_parquet_kraken2.py):
   call `kraken2-inspect` against the configured production
   Kraken2 DB and emit `all_viruses_vs_kraken2.tsv` next to the
   parquet, plus add overlap counters to `build_stats.json`.
   When `KRAKEN_DB_FOR_COMPARE` is left pointing at the
   downloaded prebuilt DB (the default), the sidecar quantifies
   the gap between the prebuilt snapshot and the newly built
   `kraken2_refseq_viral/` DB; update the key after the first
   successful local build to compare your own snapshot against
   itself as a regression check.
8. **Refresh the BLAST viral DB tarballs** via
   `update_blastdb.pl` (BLAST+ ships the script). Fetches
   `ref_viruses_rep_genomes`, `mito_rna_db` and `taxdb` from
   NCBI's pre-built tarballs into a sibling directory next to the
   parquet (`blast_refseq_viral/`), writes a `viral_rna_mito.nal`
   alias so the main pipeline's `BLASTN_DB` points at a single
   prefix querying both viral and mito/rRNA references, and
   records each tarball's fetch date in a `snapshot.tsv` manifest
   for the cross-DB-coordination audit. This rule supersedes the
   previously-manual BLAST refresh step, so all four classifier
   DBs (Kraken2, Kaiju, BLAST, taxdump-driven parquet) now share
   the same snapshot.

## Expected outputs

After a successful refresh the parquet's directory looks like:

```
/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/
├── all_viruses.parquet               (~250 MB)
├── all_viruses_build_stats.json
├── all_viruses_vs_kraken2.tsv        (one row per tax_id; ~1 MB)
├── nodes.dmp                         (~215 MB)
├── kaiju_refseq_viral/
│   ├── kaiju_refseq_viral.fmi        (~250 MB)
│   ├── nodes.dmp
│   └── names.dmp
├── kraken2_refseq_viral/
│   ├── hash.k2d                       (~600 MB)
│   ├── taxo.k2d
│   ├── opts.k2d
│   ├── seqid2taxid.map
│   └── inspect.txt                    (taxid roster sidecar)
└── blast_refseq_viral/
    ├── viral_rna_mito.nal             (the alias the main pipeline points at)
    ├── ref_viruses_rep_genomes.*      (BLAST index files from NCBI's tarball)
    ├── mito_rna_db.*                  (BLAST index files from NCBI's tarball)
    ├── taxdb.bti, taxdb.btd           (BLAST taxonomy lookup)
    └── snapshot.tsv                   (one row per tarball: name + fetch UTC)
```

`build_stats.json` records:

- `input_records` / `output_records` — input FASTA size vs after
  one-rep-per-taxid dedup
- `unique_taxids` — typically ~14 900 from viral RefSeq
- `rank_distribution` — `no rank`, `species`, `serotype` etc.
- `with_genus_taxid_count` — rows whose taxdump walk-up resolved
  to a genus (typically ~86%)
- `intersection_count` — tax_ids in **both** the parquet and the
  Kraken2 DB
- `parquet_only_count` / `kraken2_only_count` — the residual
  asymmetry. Kraken2 typically knows ~16 000 more tax_ids than the
  parquet has references for (NCBI taxonomy propagation), and the
  parquet has ~400 RefSeq-only tax_ids that the pre-built Kraken2
  tarball does not.

## Point the pipeline at the refreshed databases

Update your run config (e.g. `config/config.local.yaml`):

```yaml
VIRUS_PARQUET: "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet"
TAXDUMP_NODES: "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/nodes.dmp"
KAIJU_DB:      "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/kaiju_refseq_viral"
KRAKEN_DB:     "/Volumes/LaCie/REGIONEN/ref_dbs/KRAKEN_DB/k2_viral_20260226"
BLASTN_DB:     "/Volumes/LaCie/REGIONEN/ref_dbs/BLAST_DB/blast_db/viral_rna_mito"
```

Force `bwa_align_to_kraken_hits` and `kaiju` to re-run on the
samples you care about so the coverage references and the Kaiju
classifications pick up the new databases:

```bash
snakemake --sdm conda --cores 4 \
    --configfile config/config.local.yaml \
    --forcerun bwa_align_to_kraken_hits kaiju
```

## BLAST viral DB refresh (now automated)

The `refresh_blast` rule drives `update_blastdb.pl` against
NCBI's pre-built tarballs (`ref_viruses_rep_genomes`,
`mito_rna_db`, `taxdb`) and publishes them next to the parquet
under `blast_refseq_viral/`, with a generated
`viral_rna_mito.nal` alias the main pipeline's `BLASTN_DB`
points at. The rule also writes a `snapshot.tsv` manifest so
the operator can audit at a glance whether the BLAST refresh
co-dates with the parquet rebuild.

Override the default DB list in `refresh/config.local.yaml` if
your local install uses a different alias mix:

```yaml
BLAST_NAMES: ["ref_viruses_rep_genomes", "mito_rna_db", "taxdb"]
BLAST_PUBLISH_DIR: "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/blast_refseq_viral"
```

After the refresh, point the main pipeline at the new alias:

```yaml
BLASTN_DB: "/Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/blast_refseq_viral/viral_rna_mito"
```

## Troubleshooting

**`utime() ... No such file or directory` in a download log.**
The fix is already in: `curl` replaced `wget` for the large
download rules in this workflow. If you see it on a different
volume, the same `curl --silent --show-error --fail --location`
invocation in `refresh/refresh_virus_parquet.smk` is the pattern
to follow.

**`IncompleteFilesException` or `LockException` on retry.** Pass
`--rerun-incomplete` and (if a previous run was killed) the
workflow may also need `snakemake -s refresh/... --unlock`
once. Outputs that survived the previous attempt are kept; the
broken file gets re-downloaded.

**Kaiju classifies reads to plant / fungal taxa.** That means
the Kaiju FMI was built with the wrong header format. Verify
with `head -1 _refresh_workdir/kaiju_refseq_viral.faa` — the
header should be a bare integer (`>1980428`) not the
`>kaiju|<taxid>|<accession>` format documented on some
third-party pages. The current
[`reformat_kaiju_headers.py`](../scripts/reformat_kaiju_headers.py)
emits the bare-integer format; rerun the build if the FMI on
disk pre-dates this fix.

**Kraken2 overlap reports 0 taxids.** That means
`kraken2-inspect` was invoked with `--skip-counts`, which
suppresses the per-taxon table the parser needs. The current
[`compare_parquet_kraken2.py`](../scripts/compare_parquet_kraken2.py)
does not pass that flag; if you see a stale result rerun the
`compare_with_kraken2` rule with `--forcerun`.
