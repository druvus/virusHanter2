# Provenance: databases and application versions

Every run records which reference databases and which application
versions produced it, so a result can be traced to the exact inputs
that generated it. The record appears in three places, all derived from
the same source of truth.

## What is captured

**Application (tool) versions -- what actually ran.** The conda envs pin
only minimum versions (`kraken2>=2.1`, ...), so the resolved version is
decided when the env is materialised. A per-env probe rule reuses the
same materialised conda prefix as the real work and reads its
`conda-meta` directory, so the recorded version is exactly the one that
ran, not a lower bound. Optional stages (geNomad, QUAST, hostile,
MultiQC) are probed only when they run.

**Database build identity -- a robust stamp, not a bare mtime.** Each
reference database resolves to a build identity in preference order:

1. a `build_stats.json` sidecar written by the refresh workflow
   (`build_date_utc` + `source`), used for the parquet and, after the
   next rebuild, the Kaiju and Kraken2 viral DBs;
2. the version-bearing directory or prefix name (`checkv-db-v1.5`,
   `k2_pluspf_20240112`);
3. a representative-file mtime as the date of last resort.

Paths are never shown in full: only the parent folder and leaf are kept
(`checkv/checkv-db-v1.5`), so the operator's filesystem layout does not
leak into a shared report.

## Where it appears

- **`run_provenance_<batch>.json`** (with a flat `.tsv` companion) -- the
  machine-readable sidecar, diffable across runs. Written by
  `rules/provenance.smk::write_provenance` via
  `scripts/write_provenance.py`. This is the contract the report renders.
- **`software_versions.tsv`** -- the full resolved package table
  (`env`, `package`, `version`, `build`) merged from the probes by
  `scripts/collect_software_versions.py`.
- **`run_information_<batch>.csv`** -- trailing additive columns
  `databases_build_identity` and `tool_versions` (alongside the existing
  `databases_used`, `databases_provenance`, `databases_span_days`,
  `reporthanter_version`). These follow the 14 parity-locked columns and
  are dropped before any parity diff.
- **The per-sample HTML report** -- a `Provenance` tab (reportHanter's
  `ProvenanceSection`, fed via the `--provenance_file` CLI flag) listing
  the run scalars, the reference databases and the resolved tool
  versions.

## Cross-database coordination

`databases_span_days` reports the span between the oldest and newest DB
build dates. A span over 180 days prints a warning during aggregation:
the classifier DBs likely came from divergent snapshots, so rebuild
`VIRUS_PARQUET` and the classifier DBs from one coordinated snapshot via
`refresh/refresh_virus_parquet.smk` (see
[REFRESH_TUTORIAL.md](REFRESH_TUTORIAL.md)).
