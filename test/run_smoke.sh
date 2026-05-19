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

# Apple Silicon detection. Historically this script forced
# CONDA_SUBDIR=osx-64 because several bioconda tools (bwa, samtools, ...)
# only shipped osx-64 builds, but the Rosetta path then segfaulted on
# AVX/SIMD-heavy tools (MEGAHIT, kaiju). The upstream bioconda packages
# we depend on now ship native osx-arm64 (or noarch) builds, so we no
# longer pin the subdir and let conda choose. Kraken2 with the pluspf
# database still needs more RAM than a typical laptop has; on Apple
# Silicon prefer the standard DB or use Linux. QUAST currently lacks
# an osx-arm64 build, so flip QUAST: "FALSE" or set
# CONDA_SUBDIR=osx-64 if you enable the assembly-QC stage on a Mac.
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    echo "[smoke] Apple Silicon detected; using native osx-arm64 envs."
    echo "[smoke] If a rule fails with 'PackagesNotFoundError', that env"
    echo "[smoke] still lacks an arm64 build and may need CONDA_SUBDIR=osx-64."
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

    if [[ "$checkv_ready" == "1" ]]; then
        echo "[smoke] CheckV DB present; running full pipeline"
        snakemake --sdm conda --cores 2 --configfile "$CONFIG"
        local sample_dir="test/results/test"
        local html="$sample_dir/test_R/REPORT/test_R.html"
        local run_csv="$sample_dir/run_information_test.csv"
        [[ -s "$html"     ]] || { echo "[smoke] FAIL: missing $html";     exit 1; }
        [[ -s "$run_csv"  ]] || { echo "[smoke] FAIL: missing $run_csv";  exit 1; }
        echo "[smoke] OK: $html, $run_csv"
    else
        echo "[smoke] CheckV DB is stubbed; running everything that doesn't depend on CheckV"
        # Terminal rules whose combined dependencies cover the assembly
        # branch, the classification branch, and the coverage branch.
        # CheckV / merge_checkv_blastn / generate_report /
        # aggregate_run_information stay out of this set.
        snakemake --sdm conda --cores 2 --configfile "$CONFIG" \
            --until blastn mosdepth_kraken_hits kaiju_to_table
        local sd="test/results/test/test_R"
        local expected=(
            "$sd/BLASTN/test_R.contigs.blastn.csv"
            "$sd/KAIJU/test_R.kaiju.table.tsv"
            "$sd/KRAKEN/test_R.kraken.csv"
            "$sd/MOSDEPTH/test_R.regions.bed.gz"
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
