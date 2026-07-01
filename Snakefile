# Snakefile

# Import necessary modules
import os
from pathlib import Path

# Load configuration
configfile: "config/config.yaml"

# Validate the config against its JSON schema before any rule depends on it.
# Catches the common mistake of running the workflow against the template
# config.yaml without filling in the database paths.
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate
validate(config, "config/config.schema.yaml")


def _check_db_paths_exist(cfg: dict) -> None:
    """Fail fast with a single clean message if any required reference
    database is not on disk where the config claims.

    Catches the common operational mishap of an unmounted external
    drive: without this, the first rule that touches the missing path
    crashes mid-run with a tool-specific error (samtools, blastn,
    etc.) instead of a clear "your database is missing" report. The
    check is best-effort -- it tolerates BWA / BLAST prefix-only paths
    by checking the parent directory.
    """
    missing: list[str] = []

    def check_dir(key: str) -> None:
        path = cfg.get(key)
        if not path:
            return
        if not Path(path).is_dir():
            missing.append(f"  {key}: {path}  (expected directory)")

    def check_file(key: str) -> None:
        path = cfg.get(key)
        if not path:
            return
        if not Path(path).is_file():
            missing.append(f"  {key}: {path}  (expected file)")

    def check_prefix(key: str) -> None:
        # BWA / BLAST databases are referenced by a prefix that has no
        # file on disk; the indexed siblings live in the same dir.
        # Treat the parent dir's existence as proxy.
        path = cfg.get(key)
        if not path:
            return
        if not Path(path).parent.is_dir():
            missing.append(
                f"  {key}: {path}  (parent directory missing)"
            )

    if cfg.get("TAXDUMP_NODES"):
        check_file("TAXDUMP_NODES")
    check_prefix("HUMAN_INDEX")
    check_dir("KAIJU_DB")
    check_dir("KRAKEN_DB")
    check_prefix("BLASTN_DB")
    check_dir("CHECKV_DB")
    check_file("VIRUS_PARQUET")
    if cfg.get("GENOMAD", "FALSE") == "TRUE":
        check_dir("GENOMAD_DB")
    if cfg.get("SECONDARY_HOST_INDEX"):
        check_prefix("SECONDARY_HOST_INDEX")

    if missing:
        raise WorkflowError(
            "Reference databases are missing or unreachable. "
            "If they live on an external volume, check that it is "
            "mounted before retrying.\n" + "\n".join(missing)
        )


def _is_non_executing_invocation() -> bool:
    """True when snakemake will only build the DAG, not run any rule.

    Dry-runs (``-n`` / ``--dry-run`` / ``--dryrun``) and lint
    (``--lint``) parse the workflow without reading a single input file,
    so the reference databases need not be present on disk. The CI smoke
    job runs ``snakemake -n`` without materialising the test databases;
    gating the existence check on this keeps that dry-run green while a
    real run against an unmounted volume still fails fast. Detected from
    ``sys.argv`` because this runs at Snakefile parse time, before
    snakemake exposes a resolved settings object.
    """
    import sys

    # Snakemake's no-argument short flags that can be bundled with the
    # dry-run 'n' (printshellcmds, reason, forceall, keep-going, quiet,
    # touch). A real-run bundle that carries a value -- e.g. "-sSnakefile"
    # whose value happens to contain an 'n' -- must NOT be mistaken for a
    # dry-run, or the database existence guard would be skipped on a real
    # run against an unmounted volume.
    _BOOL_SHORT_FLAGS = set("nprFkqt")

    for arg in sys.argv:
        if arg in ("-n", "--dry-run", "--dryrun", "--lint"):
            return True
        # Treat a single-dash token as a dry-run bundle only when every
        # character is a known boolean short flag and 'n' is among them.
        if len(arg) >= 2 and arg[0] == "-" and arg[1] != "-":
            body = arg[1:]
            if "n" in body and set(body) <= _BOOL_SHORT_FLAGS:
                return True
    return False


# Skip the on-disk database check when snakemake is only building the
# DAG (dry-run / lint): those modes never touch a database, and CI runs
# the dry-run without the test fixtures present. A real run still gets
# the fail-fast guard against an unmounted reference volume.
if not _is_non_executing_invocation():
    _check_db_paths_exist(config)


def _warn_db_snapshot_mismatch(cfg: dict) -> None:
    """Emit a warning when reference databases appear to come from different
    NCBI snapshots.

    Compares build dates where a ``build_stats.json`` sidecar exists next to
    ``VIRUS_PARQUET``; otherwise falls back to directory / file modification
    times.  A spread greater than 30 days between any pair triggers the
    warning.  The check is advisory only -- the workflow is not aborted.
    """
    import json
    import sys
    import time

    MAX_SPREAD_DAYS = 30

    def _mtime(path_str: str) -> float:
        """Return the mtime of a path (file or directory) in seconds."""
        try:
            return Path(path_str).stat().st_mtime
        except OSError:
            return 0.0

    # Collect (label, epoch_seconds) for each DB.
    timestamps: list[tuple[str, float]] = []

    # VIRUS_PARQUET: prefer build_stats.json build_date_utc over mtime.
    parquet_path = cfg.get("VIRUS_PARQUET", "")
    if parquet_path:
        stats_path = Path(parquet_path).with_name(
            Path(parquet_path).stem + "_build_stats.json"
        )
        ts = 0.0
        if stats_path.is_file():
            try:
                with open(stats_path) as fh:
                    stats = json.load(fh)
                # build_date_utc is ISO-8601, e.g. "2026-05-17T14:22:33+00:00"
                from datetime import datetime, timezone
                ts = datetime.fromisoformat(
                    stats["build_date_utc"]
                ).astimezone(timezone.utc).timestamp()
            except Exception:
                ts = _mtime(parquet_path)
        else:
            ts = _mtime(parquet_path)
        if ts:
            timestamps.append(("VIRUS_PARQUET", ts))

    # KRAKEN_DB, KAIJU_DB, BLASTN_DB: directory / parent-directory mtime.
    for key in ("KRAKEN_DB", "KAIJU_DB"):
        path = cfg.get(key, "")
        if path:
            ts = _mtime(path)
            if ts:
                timestamps.append((key, ts))

    blast_db = cfg.get("BLASTN_DB", "")
    if blast_db:
        ts = _mtime(str(Path(blast_db).parent))
        if ts:
            timestamps.append(("BLASTN_DB", ts))

    if len(timestamps) < 2:
        return

    min_ts = min(t for _, t in timestamps)
    max_ts = max(t for _, t in timestamps)
    spread_days = (max_ts - min_ts) / 86400.0

    if spread_days > MAX_SPREAD_DAYS:
        oldest = min(timestamps, key=lambda x: x[1])
        newest = max(timestamps, key=lambda x: x[1])
        print(
            f"WARNING: reference databases may come from different NCBI snapshots "
            f"(spread {spread_days:.0f} days). "
            f"Oldest: {oldest[0]} "
            f"({time.strftime('%Y-%m-%d', time.gmtime(oldest[1]))}), "
            f"newest: {newest[0]} "
            f"({time.strftime('%Y-%m-%d', time.gmtime(newest[1]))}). "
            "Run refresh/refresh_virus_parquet.smk to rebuild all databases "
            "from the same snapshot.",
            file=sys.stderr,
        )


_warn_db_snapshot_mismatch(config)


# Import pipeline-side helpers. Report rendering and parsing live in the
# reportHanter package; this Snakefile only owns the data-processing side.
from scripts.functions import (
    read_file_as_blob,
    common_suffix,
    paired_reads,
    kaiju_db_files,
    fastx_file_to_df,
    wrangle_kraken,
    run_blastn,
)

# Set sample information
SAMPLES = paired_reads(config["SAMPLES"])
SUFFIX = common_suffix(config["SAMPLES"])

SAMPLES_FOLDER = config["SAMPLES"]
RESULT_FOLDER = os.path.join(config["RESULTS_FOLDER"], Path(SAMPLES_FOLDER).name)

# Determine if secondary host is specified
SECONDARY_HOST = config.get("SECONDARY_HOST_INDEX", "")
SECONDARY_HOST_NAME = config.get("SECONDARY_HOST_NAME", "")
SECONDARY_HOST_OR_NOT = bool(SECONDARY_HOST)

# Set clean list based on configuration
clean_list = [f"{RESULT_FOLDER}/analysis_done.txt"] if config.get("CLEAN", "FALSE") == "TRUE" else []

# Optional run-level QC. Default on; set MULTIQC: "FALSE" to skip.
RUN_MULTIQC = config.get("MULTIQC", "TRUE") == "TRUE"

# Optional geNomad second viral-contig classifier. Default off so the
# parity invariant holds; flip GENOMAD: "TRUE" and populate GENOMAD_DB
# to opt in.
RUN_GENOMAD_WF = config.get("GENOMAD", "FALSE") == "TRUE"

# Optional QUAST assembly assessment. Default off so the parity
# invariant holds; flip QUAST: "TRUE" to opt in. Output is also fed
# into MultiQC when both are enabled.
RUN_QUAST_WF = config.get("QUAST", "FALSE") == "TRUE"

# Active de novo assemblers. Each name in ASSEMBLERS drives one
# independent pipeline through Pilon / BLASTN / CheckV / geNomad /
# QUAST; outputs accumulate under
# `{RESULT_FOLDER}/{sample}/{assembler}/`. Default is both MEGAHIT
# and metaSPAdes. To recover byte-identical parity with the
# original virusHanter, set ASSEMBLERS: ["MEGAHIT"] in config.
ASSEMBLERS = list(
    config.get("ASSEMBLERS", ["MEGAHIT", "metaSPAdes", "rnaviralSPAdes"])
)
_VALID_ASSEMBLERS = {"MEGAHIT", "metaSPAdes", "rnaviralSPAdes"}
# Migration: the older "SPAdes" alias was renamed to "metaSPAdes"
# (it always invoked `spades.py --meta`; the new name is explicit
# about that). Catch the deprecated entry with a clear error so
# operators with old configs see exactly what to change.
if "SPAdes" in ASSEMBLERS:
    raise WorkflowError(
        "config[ASSEMBLERS] contains the deprecated alias 'SPAdes'. "
        "Rename it to 'metaSPAdes' (the rule has always run "
        "`spades.py --meta`; the new name spells that out). "
        "Existing on-disk results under {sample}/SPAdes/ are still "
        "readable; rename or copy them to {sample}/metaSPAdes/ if "
        "you want to keep them."
    )
_invalid = [a for a in ASSEMBLERS if a not in _VALID_ASSEMBLERS]
if _invalid:
    raise WorkflowError(
        "Unknown assembler(s) in config[ASSEMBLERS]: "
        + ", ".join(_invalid)
        + ". Valid choices are: "
        + ", ".join(sorted(_VALID_ASSEMBLERS))
    )

# Classifier sources that contribute taxids to the BWA reference
# set used by mosdepth coverage. Union of the three by default;
# set COVERAGE_SOURCES: ["KRAKEN"] in config to recover the
# pre-multi-source behaviour.
COVERAGE_SOURCES = list(
    config.get("COVERAGE_SOURCES", ["KRAKEN", "KAIJU", "BLAST"])
)
_VALID_COVERAGE_SOURCES = {"KRAKEN", "KAIJU", "BLAST"}
_invalid_sources = [s for s in COVERAGE_SOURCES if s not in _VALID_COVERAGE_SOURCES]
if _invalid_sources:
    raise WorkflowError(
        "Unknown coverage source(s) in config[COVERAGE_SOURCES]: "
        + ", ".join(_invalid_sources)
        + ". Valid choices are: "
        + ", ".join(sorted(_VALID_COVERAGE_SOURCES))
    )
COVERAGE_TOP_N = int(config.get("COVERAGE_TOP_N", 20))

# Optional taxdump-driven rank filter and genus walk-up. Both are
# silent no-ops when TAXDUMP_NODES is empty; the rule still runs
# with the bare multi-source union behaviour.
TAXDUMP_NODES = config.get("TAXDUMP_NODES", "") or ""
COVERAGE_RANK_FILTER = list(
    config.get(
        "COVERAGE_RANK_FILTER",
        [
            # NCBI's pseudo-root for the viral subtree. Without it,
            # the bare "Viruses" taxid (10239) propagates from
            # Kraken into the unmapped sidecar on every clinical
            # sample because viral RefSeq does not carry a
            # reference for the root.
            "acellular root",
            "realm",
            "kingdom",
            "subkingdom",
            "phylum",
            "subphylum",
            "class",
            "subclass",
            "order",
            "suborder",
            "family",
            "subfamily",
        ],
    )
)
COVERAGE_GENUS_WALKUP = config.get("COVERAGE_GENUS_WALKUP", "TRUE") == "TRUE"

# Host-removal backend: "bwa" (parity default) or "hostile" (T2T-CHM13
# via minimap2; opt-in for clinical samples where the extra
# telomeric/pericentromeric host removal matters). The pre_processing
# module owns the backend dispatch table `_HOST_BACKENDS`, validates
# `config[HOST_REMOVAL]` against it, and exposes the helpers
# `host_removed_r1` / `host_removed_r2` / `host_flagstat` that resolve
# to the active backend's outputs.

# Include rule files
include: "rules/pre_processing.smk"
include: "rules/classification.smk"
include: "rules/assembly.smk"
include: "rules/post_processing.smk"
include: "rules/provenance.smk"

# Trivial rules that should run on the submission host rather than be queued.
# Note: `wrangle_kraken`, `wrangle_pilon`, and `merge_checkv_blastn` need
# pandas/pyfastx from their `conda:` envs, and Snakemake silently ignores
# `conda:` on a localrule. They run as normal jobs.
localrules:
    all,
    clean_everything,

# Define the final targets of the workflow
rule all:
    input:
        expand(f"{RESULT_FOLDER}/{{sample}}/REPORT/{{sample}}.html", sample=SAMPLES),
        f"{RESULT_FOLDER}/run_information_{Path(SAMPLES_FOLDER).name}.csv",
        # Per-(sample, virus) detail CSV for the collaborator; concatenated
        # across samples by the aggregate_per_virus rule.
        f"{RESULT_FOLDER}/per_virus_{Path(SAMPLES_FOLDER).name}.csv",
        # Provenance: resolved tool versions + the run provenance sidecar
        # (DB build identity + versions) that the reports render. Always
        # produced so every run records what databases and applications
        # were used.
        f"{RESULT_FOLDER}/software_versions.tsv",
        f"{RESULT_FOLDER}/run_provenance_{Path(SAMPLES_FOLDER).name}.json",
        # Per-sample additive QC outputs (do not feed any other rule).
        # markdup stats exist only for the bwa backend; the hostile
        # backend has no bwa human BAM to mark duplicates on, so demanding
        # them there would force the entire redundant bwa chain to run.
        (
            expand(f"{RESULT_FOLDER}/{{sample}}/logs/human_markdup_stats.txt", sample=SAMPLES)
            if HOST_REMOVAL == "bwa"
            else []
        ),
        expand(f"{RESULT_FOLDER}/{{sample}}/MOSDEPTH/{{sample}}.mosdepth.summary.txt", sample=SAMPLES),
        # Run-level QC, gated by MULTIQC config flag (default TRUE).
        [f"{RESULT_FOLDER}/multiqc_report.html"] if RUN_MULTIQC else [],
        # Optional geNomad classifier, gated by GENOMAD config flag.
        # Per-assembler outputs accumulate under
        # `{sample}/{assembler}/GENOMAD/`.
        (
            expand(
                f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/GENOMAD/{{sample}}_improved_contigs_summary/{{sample}}_improved_contigs_virus_summary.tsv",
                sample=SAMPLES,
                assembler=ASSEMBLERS,
            )
            if RUN_GENOMAD_WF
            else []
        ),
        # Optional QUAST assembly assessment, gated by QUAST config flag.
        # One QUAST run per (sample, assembler).
        (
            expand(
                f"{RESULT_FOLDER}/{{sample}}/{{assembler}}/QUAST/report.tsv",
                sample=SAMPLES,
                assembler=ASSEMBLERS,
            )
            if RUN_QUAST_WF
            else []
        ),
        clean_list,