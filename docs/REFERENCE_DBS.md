# Reference databases

The pipeline consumes several reference databases. All currently live
under `/Volumes/LaCie/REGIONEN/ref_dbs/` and are referenced by the
production config at
`virusHanter2/config/config.production.yaml`.

| DB | Path | Used by | Refresh cadence |
|---|---|---|---|
| Human BWA index | `BWA_GENCODE_GRCH38/human_gencode*` | `bwa_human` | stable, very rarely |
| Kraken2 standard | `KRAKEN_DB/standard/` | `kraken` (lower RAM) | annual |
| Kraken2 pluspf | `KRAKEN_DB/pluspf/` | `kraken` (richer, ~82 GB hash) | annual |
| Kaiju refseq | `KAIJU_DB/refseq/` | `kaiju` | annual |
| CheckV v1.5 | `CHECKV_DB/checkv-db-v1.5/` | `checkv` | when CheckV releases a new DB |
| BLAST viral + mito | `BLAST_DB/blast_db/` + alias `viral_rna_mito.nal` | `blastn` | quarterly |
| Viral RefSeq FASTA | `VIRUS_FASTA/viral_refseq_<YYYYMMDD>.fna` | source for `all_viruses.parquet` | quarterly |
| `nucl_gb.accession2taxid.gz` | `INDIVIDUAL_VIRUS_FASTA/nucl_gb.accession2taxid.gz` | source for `all_viruses.parquet` | with each viral RefSeq refresh |
| `all_viruses.parquet` | `INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet` | `bwa_align_to_kraken_hits`, `per_virus_metrics` | rebuild after FASTA refresh |
| geNomad DB (optional) | `GENOMAD_DB/genomad_db/` | `genomad` (only when `GENOMAD: "TRUE"`) | when geNomad releases a new DB |

## CheckV on a Mac-mounted external volume

If the CheckV DB is hosted on an external drive that has ever been
written to by macOS Finder (HFS+/APFS volumes formatted by macOS, or
network shares browsed from macOS), `hmm_db/checkv_hmms/` will end up
with AppleDouble companion files (`._1.hmm`, `._2.hmm`, ...). These
are tiny resource-fork metadata stubs that macOS creates alongside
real files on volumes that lack native extended-attribute support.

CheckV's hmmsearch driver lists the HMM directory and feeds every
entry to its multiprocessing pool. The ghost files are not valid
HMMs, so hmmsearch errors on each one with
`File existence/permissions problem in trying to open HMM file ...`,
and CheckV reports a generic
`Error: 80 hmmsearch tasks failed. Program should be rerun.` even
when every real HMM ran cleanly.

Strip them once after each refresh:

```
find /path/to/checkv-db-v1.5 -name '._*' -delete
```

For a permanent fix, host the CheckV DB on an APFS-local volume so
Finder never creates the shadows in the first place.

## Rebuilding `all_viruses.parquet`

The parquet is keyed by `(name, sequence, tax_id)` and is consumed by
`bwa_align_to_kraken_hits` to pick reference sequences for the top
Kraken2 viral taxids, and by `per_virus_metrics` to attribute
contigs and aggregate coverage per taxid. Refresh it whenever the
viral RefSeq FASTA changes, using the bounded-memory builder.

```
cd /Users/andreassjodin/Code/regionen/virusHanter2
python scripts/build_virus_parquet.py \
    --fasta /Volumes/LaCie/REGIONEN/ref_dbs/VIRUS_FASTA/viral_refseq_<YYYYMMDD>.fna \
    --taxid /Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/nucl_gb.accession2taxid.gz \
    --out   /Volumes/LaCie/REGIONEN/ref_dbs/INDIVIDUAL_VIRUS_FASTA/all_viruses.parquet
```

The script streams the gzipped `nucl_gb.accession2taxid` and keeps
only the entries whose accession appears in the viral FASTA (~30k
keys). Memory stays bounded.

Last successful rebuild on this workstation (2026-05-17): 19,149
records, 14,899 unique tax_ids, 39 records (0.2%) without a taxid,
mean reference length 30 kb. The 39 unmatched records are normal —
typically very recent submissions whose accession2taxid mapping has
not propagated yet.

## Sources for the inputs

- Viral RefSeq FASTA: download from
  `https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/` and concatenate
  the `*.genomic.fna.gz` files into one FASTA named
  `viral_refseq_<YYYYMMDD>.fna`.
- `nucl_gb.accession2taxid.gz`: download from
  `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/`
  along with its `.md5` for verification.
- Kraken2 indices: download from
  `https://genome-idx.s3.amazonaws.com/kraken/` and extract in
  place. The pluspf 2024-01-12 tarball is ~63 GB compressed and
  expands to ~96 GB.
- Kaiju indices: download from `https://kaiju.binf.ku.dk/server` or
  build with `kaiju-makedb`. The refseq build is ~22 GB.
- CheckV v1.5: `https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz`.
- BLAST viral DBs: refresh `ref_viruses_rep_genomes`, `mito_rna_db`,
  and the `taxdb` files via `update_blastdb.pl` from the BLAST+
  toolkit. The `viral_rna_mito.nal` alias is hand-written and
  references both the viral and mito BLAST databases.
- geNomad DB (only if you opt into the genomad rule): the
  recommended fetch is the geNomad CLI itself, which fetches the
  current release into a `genomad_db/` subdirectory of the path
  you give it:

  ```bash
  conda activate <env with genomad>
  cd /Volumes/LaCie/REGIONEN/ref_dbs
  mkdir -p GENOMAD_DB && cd GENOMAD_DB
  genomad download-database .
  # produces ./genomad_db/
  ```

  Then set `GENOMAD_DB: "/Volumes/LaCie/REGIONEN/ref_dbs/GENOMAD_DB/genomad_db"`
  in your config. A manual fallback is the Zenodo mirror at
  `https://zenodo.org/records/14886553`; extract the tarball and
  point `GENOMAD_DB` at the resulting `genomad_db/` directory.

  Performance caveat: `genomad end-to-end` calls mmseqs2, which
  does many small random reads against the ~228 k marker profiles.
  Running with the DB on an external USB drive is impractical
  (9+ hour stall observed on this workstation with the geNomad DB
  on LaCie). For real-data runs, stage the DB on local SSD or a
  fast network filesystem before flipping `GENOMAD: "TRUE"`.

## Apple Silicon / RAM-limited host notes

- Kraken2 with the pluspf hash needs more RAM than is available on a
  typical laptop. Use the `standard` build for local debugging; run
  pluspf on a Linux host with >= 96 GB RAM.
- Kaiju's refseq `.fmi` is ~22 GB; it loads entirely into memory and
  the same RAM constraint applies.
- The viral parquet, BWA indices, and CheckV / mosdepth / fastp all
  run within ~16 GB.
