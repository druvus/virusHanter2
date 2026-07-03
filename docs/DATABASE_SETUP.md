# Database setup

This document is the single entry point for obtaining and configuring
all reference databases consumed by `virusHanter2`. Read it before
your first run and again whenever you refresh the classification stack.

For background on why snapshot co-ordination matters, see
[REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md). For the full list of
config keys, see [CONFIGURATION.md](CONFIGURATION.md).
For per-database disk paths used in production, see
[REFERENCE_DBS.md](REFERENCE_DBS.md).

> **Layout note.** The `/path/to/<TOOL>_DB/...` examples here describe
> obtaining each database individually (the production per-tool layout).
> For a fresh viral-only deploy, the refresh workflow builds the viral
> classifier DBs together under `$VH2_ROOT/refdbs/virus_ref/`
> (`all_viruses.parquet`, `nodes.dmp`, `kaiju_refseq_viral/`,
> `kraken2_refseq_viral/`, `blast_refseq_viral/`), with `human/` and
> `checkv/` alongside -- see [DEPLOY_LINUX.md](DEPLOY_LINUX.md) and
> [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md). Use one layout consistently.

---

## Resource requirements at a glance

Sizes are approximate. Compressed download sizes are given where a
tarball is the normal delivery mechanism; expanded (on-disk) sizes
are what the pipeline actually sees at runtime. RAM figures are the
working-set peak for the step that reads the index into memory.

| Database          | Config key        | Required?                  | Disk (expanded) | RAM peak  | Download  |
|-------------------|-------------------|----------------------------|-----------------|-----------|-----------|
| Human BWA index   | `HUMAN_INDEX`     | Yes (bwa host removal)     | ~14 GB          | ~4 GB     | ~3.2 GB   |
| Kraken2 standard  | `KRAKEN_DB`       | Yes (one of the two)       | ~30 GB          | ~30 GB    | ~20 GB    |
| Kraken2 pluspf    | `KRAKEN_DB`       | Yes (one of the two)       | ~96 GB          | ~96 GB    | ~63 GB    |
| Kaiju refseq      | `KAIJU_DB`        | Yes (one of the two)       | ~22 GB          | ~24 GB    | ~22 GB    |
| Kaiju viruses     | `KAIJU_DB`        | Yes (one of the two)       | ~318 MB         | ~1 GB     | ~318 MB   |
| BLAST viral alias | `BLASTN_DB`       | Yes                        | ~3 GB           | low       | ~2 GB     |
| CheckV v1.5       | `CHECKV_DB`       | Yes                        | ~2 GB           | low       | ~1.5 GB   |
| Viral parquet     | `VIRUS_PARQUET`   | Yes                        | ~250 MB         | ~1 GB     | built     |
| NCBI nodes.dmp    | `TAXDUMP_NODES`   | Strongly recommended       | ~215 MB         | low       | ~60 MB*   |
| geNomad DB        | `GENOMAD_DB`      | Only when `GENOMAD: "TRUE"`| ~4 GB           | ~6 GB     | ~4 GB     |
| hostile index     | `HOSTILE_INDEX`   | Only when `HOST_REMOVAL: "hostile"` | ~2 GB | low    | ~2 GB     |

\* `nodes.dmp` is extracted from `taxdump.tar.gz` (~60 MB compressed);
the refresh workflow handles this automatically.

**Laptop / Apple Silicon note.**
Kraken2 pluspf (~96 GB RAM) and Kaiju refseq (~22 GB RAM) both exceed
what is available on a 16--18 GB MacBook. Use Kraken2 standard and
Kaiju viruses for local development; run pluspf and refseq on a Linux
host with >= 96 GB RAM for production (see the decision notes below).

---

## Decision notes

### Kraken2: standard vs pluspf

Both databases are downloaded pre-built from
`https://genome-idx.s3.amazonaws.com/kraken/`.

- **`standard`** (~30 GB) covers bacteria, archaea, viruses and human.
  Suitable for local development and samples where bacteria/archaea
  context is not the primary interest.
- **`pluspf`** (~96 GB) adds protozoa and fungi on top of standard.
  Recommended for production clinical batches where those kingdoms may
  generate reads, and where calling out contamination from common
  environmental fungi matters. The 2024-01-12 pluspf tarball is the
  current production choice.

Set `KRAKEN_DB` in your config to the directory of whichever build you
choose. There is no pipeline-level flag to switch; the path is the
decision.

### Kaiju: refseq vs viruses

Both databases can be downloaded pre-built from
`https://kaiju.binf.ku.dk/server` or built locally with `kaiju-makedb`.

- **`refseq`** (~22 GB FMI). Covers the full NCBI RefSeq protein set.
  Recommended for all production clinical workloads. Requires a host
  with at least 24--32 GB of RAM; on a 16--18 GB laptop the index
  pages to swap and classification takes hours per sample (or OOM-kills).
- **`viruses`** (~318 MB FMI). Covers RefSeq complete viral genomes only.
  Fits comfortably on a laptop. However, content gaps make it unsuitable
  for clinical work: a 2026-05-21 audit found 0 Anelloviridae proteins
  (Torque teno virus and relatives) and only 1 Pegiviridae protein (GB
  virus C and relatives). Batches dominated by those families will
  report `kaiju_virus_percent = 0` regardless of Kraken2 signal. Use
  this database only for local development and smoke tests, and treat a
  blank Kaiju panel as expected. Additionally, the pre-built `viruses/`
  directory lacks `nodes.dmp` and `names.dmp` alongside the FMI; you
  must create symlinks before running (see gotchas below).

The alternative to the pre-built databases is the refresh workflow
(`refresh/refresh_virus_parquet.smk`), which builds a Kaiju FMI from
the same NCBI viral RefSeq snapshot as the parquet and closes the
Anelloviridae / Pegiviridae gap. See [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md).

---

## Database-by-database setup

### 1. Human BWA index (`HUMAN_INDEX`)

**Purpose.** The `bwa_human` rule aligns every sample's reads to the
human reference genome and discards mapped pairs as host. The index
prefix is passed directly to `bwa mem -k 26`.

**Config key.**
```yaml
HUMAN_INDEX: "/path/to/BWA_GENCODE_GRCH38/human_gencode"
```
The value is the BWA index *prefix* (the stem shared by the `.amb`,
`.ann`, `.bwt`, `.pac`, `.sa` files), not the FASTA path.

**Build command.** There is no pre-built BWA index for GRCh38 available
from a canonical source; you must index your own copy of the genome.
Download the GRCh38 primary assembly FASTA from GENCODE:

```bash
# Download the current GENCODE primary assembly (adjust the release
# number to the current release; check https://www.gencodegenes.org/human/).
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# Decompress (bwa index can read gzip but decompressing first is faster
# and avoids a second pass through the gzip stream on subsequent runs).
gunzip GRCh38.primary_assembly.genome.fa.gz

# Build the BWA index (wall time ~1.5 h on a single core; ~14 GB on disk).
bwa index -p human_gencode GRCh38.primary_assembly.genome.fa
```

Set `HUMAN_INDEX` to the absolute path of the `human_gencode` prefix
you chose. The FASTA itself is not read by the pipeline after indexing
and does not need to remain on the same volume.

**Gotchas.**
- The index must be built once per BWA version; BWA 0.7.17 and 0.7.18
  produce compatible indices, but mixing major releases can silently
  produce wrong alignments.
- If you have an existing GRCh38 BWA index from another tool or
  pipeline, confirm the prefix resolves all five index files before
  pointing the config at it.
- Only consulted when `HOST_REMOVAL: "bwa"` (the parity default). Not
  needed when `HOST_REMOVAL: "hostile"`.

---

### 2. Kraken2 database (`KRAKEN_DB`)

**Purpose.** Kraken2 performs k-mer-based taxonomic classification of
the host-cleaned reads. The top-N viral hits from Kraken2 (and
optionally Kaiju and BLAST) drive the reference-selection step for the
mosdepth coverage panel.

**Config key.**
```yaml
KRAKEN_DB: "/path/to/KRAKEN_DB/standard/"
# or
KRAKEN_DB: "/path/to/KRAKEN_DB/pluspf/"
```
The directory must contain `hash.k2d`, `opts.k2d`, and `taxo.k2d`.

**Source.** Pre-built archives are available at:
```
https://genome-idx.s3.amazonaws.com/kraken/
```

Download and extract in place:
```bash
# Standard build (adjust the snapshot date to the version you want).
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20231009.tar.gz
mkdir -p KRAKEN_DB/standard && tar -xzf k2_standard_20231009.tar.gz -C KRAKEN_DB/standard

# OR pluspf (larger; see decision note above).
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240112.tar.gz
mkdir -p KRAKEN_DB/pluspf && tar -xzf k2_pluspf_20240112.tar.gz -C KRAKEN_DB/pluspf
```

Alternatively, for a tightly co-ordinated snapshot (parquet + Kaiju +
Kraken2 viral + taxdump all from the same RefSeq release), use the
refresh workflow, which builds a `kraken2_refseq_viral/` index from
scratch:
```
refresh/refresh_virus_parquet.smk
```
See [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md).

**Gotchas.**
- The pluspf hash alone requires ~96 GB of free RAM; ensure the
  production host meets this before using it.
- Publicly hosted `k2_viral_*` snapshots occasionally omit individual
  genomes (e.g. the Feb 2026 snapshot lacked HSV-2 / NC_001798). The
  refresh workflow builds Kraken2 from the same RefSeq pull as the
  parquet to close this gap.

---

### 3. Kaiju database (`KAIJU_DB`)

**Purpose.** Kaiju performs protein-level classification of the
host-cleaned reads by translating them in all six frames and querying
against a Burrows-Wheeler FM-index of reference proteins. Its output
feeds the same coverage-reference selection pipeline as Kraken2.

**Config key.**
```yaml
KAIJU_DB: "/path/to/KAIJU_DB/refseq/"
# or (development / smoke only)
KAIJU_DB: "/path/to/KAIJU_DB/viruses/"
```
The directory must contain an `.fmi` file, `nodes.dmp`, and `names.dmp`.

**Source.** Pre-built databases can be downloaded from:
```
https://kaiju.binf.ku.dk/server
```
or built locally with `kaiju-makedb` from the kaiju bioconda package.

For a co-ordinated snapshot, the refresh workflow builds
`kaiju_refseq_viral.fmi` from the same viral RefSeq protein FASTA as
the parquet and closes the Anelloviridae / Pegiviridae gap in the
pre-built `viruses/` database:
```
refresh/refresh_virus_parquet.smk
```

**Gotchas.**

*Missing taxonomy files in the `viruses/` database.*
The pre-built `viruses/` directory does not include `nodes.dmp` or
`names.dmp` alongside the FMI. The pipeline expects all three files in
the same directory (the config key is a directory path, not a file
prefix). Create symlinks before running:

```bash
# Assuming the parent KAIJU_DB/ directory holds the full RefSeq
# nodes.dmp and names.dmp (true for the pre-built refseq/ layout).
cd /path/to/KAIJU_DB/viruses/
ln -s ../nodes.dmp nodes.dmp
ln -s ../names.dmp names.dmp
```

*Kaiju header format.*
If Kaiju classifies reads to plant or fungal taxa unexpectedly, the FMI
was built with the wrong header format. The header in the protein FASTA
should be a bare integer (`>1980428`), not `>kaiju|<taxid>|<accession>`.
The refresh workflow's `reformat_kaiju_headers.py` emits the correct
format; rebuild the FMI if it pre-dates that fix.

---

### 4. BLAST viral database (`BLASTN_DB`)

**Purpose.** The `blastn` rule searches assembled contigs (after Pilon
polishing) against a nucleotide BLAST database combining
`ref_viruses_rep_genomes` and `mito_rna_db`. The hits feed contig
annotation and the coverage-reference selection pipeline.

**Config key.**
```yaml
BLASTN_DB: "/path/to/BLAST_DB/blast_db/viral_rna_mito"
```
The value is the `.nal` alias *prefix* (no extension), not a directory.

**Source.** The refresh workflow (`refresh/refresh_virus_parquet.smk`)
drives `update_blastdb.pl` (shipped with the BLAST+ toolkit) to fetch
`ref_viruses_rep_genomes`, `mito_rna_db`, and `taxdb` from NCBI's
pre-built tarballs into a `blast_refseq_viral/` directory next to the
parquet and writes the `viral_rna_mito.nal` alias automatically. Use
the refresh workflow as the recommended path.

For a manual install (if you do not want to run the refresh workflow):

```bash
mkdir -p BLAST_DB/blast_db && cd BLAST_DB/blast_db

# Fetch the pre-built tarballs from NCBI (requires BLAST+ on PATH).
update_blastdb.pl --decompress ref_viruses_rep_genomes
update_blastdb.pl --decompress mito_rna_db
update_blastdb.pl --decompress taxdb

# Write a multi-volume alias that queries both databases.
cat > viral_rna_mito.nal <<'EOF'
TITLE viral_rna_mito
DBLIST ref_viruses_rep_genomes mito_rna_db
EOF
```

Then set `BLASTN_DB` to the absolute path of the `viral_rna_mito`
prefix (without `.nal`).

**Gotchas.**
- The `taxdb.bti` / `taxdb.btd` files must reside in the same
  directory as the alias; they provide the taxonomy lookup that
  annotates BLAST hits with taxids and lineages.
- Quarterly refresh cadence is recommended; the viral representative
  genomes set changes materially between NCBI releases.

---

### 5. CheckV database (`CHECKV_DB`)

**Purpose.** CheckV evaluates the completeness and contamination of
assembled viral contigs. The `checkv` rule runs `checkv end_to_end`
against this database per sample per assembler.

**Config key.**
```yaml
CHECKV_DB: "/path/to/CHECKV_DB/checkv-db-v1.5"
```

**Source.** Prefer `checkv download_database`, which fetches the
database version matching the installed CheckV (currently
`checkv-db-v1.5`):
```bash
mamba create -n checkv -c conda-forge -c bioconda 'checkv>=1.1.1'
conda run -n checkv checkv download_database /path/to/CHECKV_DB
```

Or download the tarball directly:
```bash
curl -LO https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
tar -xzf checkv-db-v1.5.tar.gz
```

The extracted directory (`checkv-db-v1.5/`) is the value of `CHECKV_DB`.

**Gotchas.**

*AppleDouble companion files (macOS volumes).*
If the CheckV database directory has ever been written to or browsed
from macOS Finder (HFS+ / APFS volumes, or network shares accessed via
Finder), the `hmm_db/checkv_hmms/` subdirectory will accumulate
AppleDouble companion files (`._1.hmm`, `._2.hmm`, ...). These are
invisible metadata stubs that macOS creates on volumes without native
extended-attribute support. CheckV's hmmsearch driver tries to open
every file in that directory, fails on each stub, and exits with
`Error: 80 hmmsearch tasks failed. Program should be rerun.` even
when every real HMM ran cleanly.

Strip the stubs once after each refresh or after any macOS access:

```bash
find /path/to/checkv-db-v1.5 -name '._*' -delete
```

For a permanent fix, keep the CheckV database on an APFS-local volume
so Finder never creates the stubs.

---

### 6. Viral reference parquet (`VIRUS_PARQUET`)

**Purpose.** A Parquet table with columns `(name, sequence, tax_id,
rank, genus_taxid)` holding one representative nucleotide sequence per
viral tax_id. The `bwa_align_to_kraken_hits` rule queries the parquet
to select reference sequences for the mosdepth coverage step, and the
`per_virus_metrics` rule uses it to attribute contigs.

**Config key.**
```yaml
VIRUS_PARQUET: "/path/to/INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet"
```

**Build command.** The parquet is not available as a download; it is
built by the refresh workflow from NCBI viral RefSeq data:

```bash
conda activate virushanter
cd virusHanter2

# First time only: copy the example config and edit the output paths.
cp refresh/config.yaml refresh/config.local.yaml
$EDITOR refresh/config.local.yaml

# Build. Re-running is idempotent; pass --rerun-incomplete after a
# partial run.
snakemake -s refresh/refresh_virus_parquet.smk \
    --configfile refresh/config.local.yaml --cores 4 \
    --sdm conda
```

The workflow downloads all `viral.*.1.genomic.fna.gz` parts from
`https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/`, deduplicates to
one representative sequence per tax_id, and enriches each row with
`rank` and `genus_taxid` from the co-downloaded taxdump. It also
publishes `nodes.dmp` next to the parquet (see `TAXDUMP_NODES` below).

For a manual legacy-mode build (viral RefSeq only, no taxdump
enrichment):

```bash
python scripts/build_virus_parquet.py \
    --source refseq --no-one-rep-per-taxid \
    --fasta /path/to/viral_refseq_<YYYYMMDD>.fna \
    --taxid /path/to/nucl_gb.accession2taxid.gz \
    --out   /path/to/all_viruses.parquet
```

Source URLs for manual downloads:
- Viral RefSeq FASTA parts:
  `https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/`
  (concatenate all `viral.*.genomic.fna.gz` parts into one FASTA)
- Accession-to-taxid mapping:
  `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz`

**Gotchas.**
- The `rank` and `genus_taxid` columns are only present when the
  parquet was built with `--taxdump-nodes nodes.dmp`. Older three-column
  parquets still work, but the rank filter and genus walk-up in
  `bwa_align_to_kraken_hits` degrade to no-ops.
- Rebuild the parquet after every RefSeq snapshot update; point the
  config at `VIRUS_PARQUET` and `TAXDUMP_NODES` simultaneously so the
  taxonomy metadata matches the sequences.

---

### 7. NCBI taxonomy nodes (`TAXDUMP_NODES`)

**Purpose.** An uncompressed NCBI `nodes.dmp` that enables three
downstream behaviours in `bwa_align_to_kraken_hits`:
- **Rank filter** (`COVERAGE_RANK_FILTER`): drops high-rank propagation
  hits (Viruses, family-level rows, etc.) that have no per-taxid
  reference sequence and would otherwise flood `unmapped_taxids.tsv`.
- **Genus walk-up** (`COVERAGE_GENUS_WALKUP`): substitutes a genus
  representative when a classifier hit is absent from the parquet.
- **ICTV species canonicalisation**: rewrites classifier and BLAST
  outputs to the current ICTV binomial species name and adds an
  `aliases` column with legacy NCBI names, acronyms, and common names.

This key is technically optional (empty string disables the above
features) but is strongly recommended for all production runs.

**Config key.**
```yaml
TAXDUMP_NODES: "/path/to/INDIVIDUAL_VIRUS_FASTA/nodes.dmp"
```

The sibling `names.dmp` must reside in the same directory; the
pipeline locates it via the parent path of `TAXDUMP_NODES`.

**Source.** The refresh workflow downloads and extracts `nodes.dmp` and
`names.dmp` automatically and publishes them next to the parquet. If
you need to extract them manually:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extract only the two files needed (avoids unpacking the full ~1 GB
# tarball contents).
tar -xzf taxdump.tar.gz nodes.dmp names.dmp
```

Set `TAXDUMP_NODES` to the absolute path of the extracted `nodes.dmp`.
Ensure `names.dmp` is in the same directory.

**Gotchas.**
- If `TAXDUMP_NODES` is set but the file does not exist, the pipeline
  fails at start-up (schema validation catches missing paths).
- Refresh the taxdump alongside the parquet and Kaiju/Kraken2 databases
  to keep the taxonomy snapshot consistent. The refresh workflow handles
  this automatically.

---

### 8. geNomad database (`GENOMAD_DB`)

**Purpose.** geNomad predicts whether assembled contigs are of viral or
plasmid origin and annotates them with marker scores. The `genomad`
rule (`GENOMAD: "TRUE"`) runs `genomad end-to-end` per sample per
assembler; its per-contig scores are appended as additive columns in
`per_virus_metrics.csv`.

This database is only required when `GENOMAD: "TRUE"`.

**Config key.**
```yaml
GENOMAD_DB: "/path/to/GENOMAD_DB/genomad_db"
```
The value is the `genomad_db/` subdirectory, not its parent.

**Source.** The recommended download uses the geNomad CLI (in any env
that has geNomad installed):

```bash
conda activate <env_with_genomad>
mkdir -p /path/to/GENOMAD_DB
cd /path/to/GENOMAD_DB
genomad download-database .
# produces ./genomad_db/
```

A manual fallback is the Zenodo mirror:
```
https://zenodo.org/records/14886553
```
Download and extract the tarball, then point `GENOMAD_DB` at the
resulting `genomad_db/` directory.

**Gotchas.**
- `genomad end-to-end` calls mmseqs2, which performs many small random
  reads against ~228 k marker profiles. Running with the database on an
  external USB drive is impractical (a 9-hour stall was observed with
  the geNomad database on an external USB drive). Stage the database on
  local SSD or a fast network filesystem before
  enabling `GENOMAD: "TRUE"`.
- The default `GENOMAD_SPLITS: 4` partitions the mmseqs search to keep
  peak memory under ~6 GB. Raise this value on memory-constrained hosts
  or lower it on hosts with abundant RAM to restore auto-split
  behaviour.

---

### 9. hostile host-removal index (`HOSTILE_INDEX`)

**Purpose.** When `HOST_REMOVAL: "hostile"`, the pipeline uses the
hostile tool (Bede et al.) instead of BWA for human-read removal.
hostile aligns reads with minimap2 against a T2T-CHM13 + HLA +
decoy reference with viral and phage regions masked, so reads from
endogenous retroviruses, phage-related host elements and other
host-embedded viral sequences survive the filter (the unmasked
`human-t2t-hla` index would drop them as host, which is incorrect for
a viral metagenomics pipeline).

This database is only required when `HOST_REMOVAL: "hostile"`.

**Config key.**
```yaml
HOSTILE_INDEX: "/path/to/HOSTILE/human-t2t-hla.rs-viral-202401_ml-phage-202401"
# Or leave empty to let hostile manage the cache itself.
HOSTILE_INDEX: ""
```

When `HOSTILE_INDEX` is empty, hostile downloads the default index
(`human-t2t-hla.rs-viral-202401_ml-phage-202401`) into its own managed
cache on first use. To pre-download to a specific location:

```bash
conda activate virushanter   # or any env with hostile installed
hostile index fetch \
    --name human-t2t-hla.rs-viral-202401_ml-phage-202401 \
    --out /path/to/HOSTILE/
```

The bundle is approximately 2 GB.

**Gotchas.**
- Use the masked variant `human-t2t-hla.rs-viral-202401_ml-phage-202401`,
  not the bare `human-t2t-hla`. The bare index lacks the viral / phage
  mask, so it treats endogenous-retroviral and phage-homologous reads
  as host and removes them — false-positive host removal that discards
  signal a viral metagenomics pipeline must retain.
- `HUMAN_INDEX` (BWA) is not consulted when `HOST_REMOVAL: "hostile"`;
  you do not need to build the BWA index if hostile is your only
  host-removal backend.

---

## Automated co-ordinated refresh

The classification databases (Kraken2 viral, Kaiju FMI, VIRUS_PARQUET,
BLAST viral alias, taxdump) should ideally be built from the same NCBI
RefSeq snapshot so that every tax_id in the parquet has a
representative in both the nucleotide and protein classifiers, and the
taxdump walk-up uses the same taxonomy version.

The `refresh/refresh_virus_parquet.smk` workflow automates this
end-to-end: one `snakemake` invocation downloads, builds, and publishes
all four databases plus the taxdump. After a successful refresh, update
your run config to point the five keys at the newly published paths:

```yaml
VIRUS_PARQUET: "/path/to/INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet"
TAXDUMP_NODES: "/path/to/INDIVIDUAL_VIRUS_FASTA/nodes.dmp"
KAIJU_DB:      "/path/to/INDIVIDUAL_VIRUS_FASTA/kaiju_refseq_viral"
KRAKEN_DB:     "/path/to/kraken2_refseq_viral"
BLASTN_DB:     "/path/to/blast_refseq_viral/viral_rna_mito"
```

See [REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md) for the full operator
workflow, expected outputs, troubleshooting, and the
`all_viruses_vs_kraken2.tsv` overlap sidecar that quantifies any
residual tax_id asymmetry between the parquet and Kraken2.

The external databases (HUMAN_INDEX, CheckV, geNomad, hostile) follow
their own release cycles and are updated independently.
