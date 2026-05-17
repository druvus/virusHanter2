#!/usr/bin/env bash
#
# Smoke runner for virusHanter2.
#
#   ./test/run_smoke.sh                # lint + dry-run only (no fixtures)
#   ./test/run_smoke.sh --build        # build mini-DBs (calls build_fixtures.sh)
#   ./test/run_smoke.sh --full         # build (if needed) + run the pipeline
#                                      # against the mini fixtures
#
# --full degrades automatically when the CheckV database has not been
# materialised (test/mini_db/checkv contains only the .stub sentinel):
# the smoke runs `snakemake --until blastn` and asserts that the BLASTN
# output exists. Provide a real CheckV DB under test/mini_db/checkv to
# extend the smoke through merge_checkv_blastn, generate_report, and
# aggregate_run_information.

set -euo pipefail

cd "$(dirname "$0")/.."

CONFIG="test/config.test.yaml"
MINI="test/mini_db"

# Many bioconda packages (bwa, samtools, etc.) do not ship an osx-arm64
# build. On Apple Silicon, fall back to osx-64 binaries via Rosetta by
# pinning CONDA_SUBDIR before any env is materialised. No-op on Linux.
#
# Known Apple-Silicon-only failures observed against real input
# (2026-05-17, sub-sampled 251015 batch, MacBook Pro):
#   - bam2plot: polars >= 1.40 segfaults via Rosetta; mitigated with
#     polars-lts-cpu, but the env still needs care.
#   - kaiju:     SIGSEGV when loading the ~22 GB refseq FMI under Rosetta.
#   - kraken2:   the standard hash.k2d (~21 GB) does not fit in RAM on
#                this host; exits with "Error reading in hash table".
#   - megahit:   megahit_core_popcnt uses popcnt/AVX instructions that
#                Rosetta does not emulate cleanly; SIGSEGV on real inputs.
# The smoke fixtures stay green because they are small enough that the
# AVX path is not exercised. Production runs must happen on Linux.
APPLE_SILICON=0
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    export CONDA_SUBDIR="${CONDA_SUBDIR:-osx-64}"
    APPLE_SILICON=1
    echo "[smoke] Apple Silicon detected; using CONDA_SUBDIR=$CONDA_SUBDIR"
    echo "[smoke] (bam2plot/polars >= 1.40 segfaults under Rosetta;"
    echo "[smoke]  the --full path will stop at merge_checkv_blastn on this host.)"
fi

build_if_requested() {
    # build_fixtures is itself a Snakemake workflow, so rerunning it is
    # cheap when outputs already exist. Always call it so missing FASTQ
    # or DB files get regenerated.
    echo "[smoke] running fixture build (Snakemake decides what to rebuild)"
    bash test/build_fixtures.sh
}

run_lint_and_dryrun() {
    # `snakemake --lint` exits non-zero whenever any style warning fires; the
    # warnings here are all style (no log: on trivial rules, run: blocks that
    # could be scripts) and do not block execution. Treat lint output as
    # advisory.
    echo "[smoke] snakemake --lint --configfile $CONFIG (advisory)"
    snakemake --lint --configfile "$CONFIG" || true
    echo "[smoke] snakemake -n --sdm conda --configfile $CONFIG"
    snakemake -n --sdm conda --configfile "$CONFIG"
}

run_full() {
    # A real CheckV database has a `genome_db` subdir. Anything else (a
    # missing directory, or one populated only with the smoke stub) is
    # treated as absent.
    local checkv_ready=0
    if [[ -d "$MINI/checkv/genome_db" ]]; then
        checkv_ready=1
    fi

    if [[ "$checkv_ready" == "1" && "$APPLE_SILICON" == "0" ]]; then
        echo "[smoke] CheckV DB present; running full pipeline"
        snakemake --sdm conda --cores 2 --configfile "$CONFIG"
        local sample_dir="test/results/test"
        local html="$sample_dir/test_R/REPORT/test_R.html"
        local run_csv="$sample_dir/run_information_test.csv"
        [[ -s "$html"     ]] || { echo "[smoke] FAIL: missing $html";     exit 1; }
        [[ -s "$run_csv"  ]] || { echo "[smoke] FAIL: missing $run_csv";  exit 1; }
        echo "[smoke] OK: $html, $run_csv"
    elif [[ "$checkv_ready" == "1" && "$APPLE_SILICON" == "1" ]]; then
        echo "[smoke] CheckV DB present, but Apple Silicon / Rosetta blocks bam2plot."
        echo "[smoke] Running everything that does not need bam2plot."
        snakemake --sdm conda --cores 2 --configfile "$CONFIG" \
            --until merge_checkv_blastn kaiju_to_table
        local merged="test/results/test/test_R/CHECKV/test_R.merged.csv"
        [[ -s "$merged" ]] || { echo "[smoke] FAIL: missing $merged"; exit 1; }
        echo "[smoke] OK: $merged"
        echo "[smoke] (bam2plot + generate_report + aggregate skipped; run on Linux for the HTML report.)"
    else
        echo "[smoke] CheckV DB is stubbed; running everything that doesn't depend on CheckV"
        # Three terminal rules whose combined dependencies cover the
        # assembly branch, the classification branch, and the coverage
        # branch. CheckV / merge_checkv_blastn / generate_report /
        # aggregate_run_information stay out of this set.
        snakemake --sdm conda --cores 2 --configfile "$CONFIG" \
            --until blastn bam2plot kaiju_to_table
        local sd="test/results/test/test_R"
        local expected=(
            "$sd/BLASTN/test_R.contigs.blastn.csv"
            "$sd/KAIJU/test_R.kaiju.table.tsv"
            "$sd/KRAKEN/test_R.kraken.csv"
            "$sd/COVERAGE_PLOTS"
        )
        for f in "${expected[@]}"; do
            [[ -e "$f" ]] || { echo "[smoke] FAIL: missing $f"; exit 1; }
        done
        echo "[smoke] OK:"
        for f in "${expected[@]}"; do echo "  - $f"; done
        echo "[smoke] (HTML report stage skipped; provide a real CheckV DB to enable it.)"
    fi
}

case "${1:-}" in
    --build)
        build_if_requested 1
        ;;
    --full)
        build_if_requested 0
        run_lint_and_dryrun
        run_full
        ;;
    "")
        run_lint_and_dryrun
        ;;
    *)
        echo "Usage: $0 [--build | --full]"
        exit 64
        ;;
esac
