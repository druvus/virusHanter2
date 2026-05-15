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

build_if_requested() {
    if [[ "$1" == "1" ]] || [[ ! -d "$MINI/human" ]]; then
        echo "[smoke] (re)building mini fixtures"
        bash test/build_fixtures.sh
    fi
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
    local checkv_ready=0
    if [[ -d "$MINI/checkv" ]] && [[ ! -e "$MINI/checkv/.stub" ]]; then
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
        echo "[smoke] CheckV DB is stubbed; running up to BLASTN only"
        snakemake --sdm conda --cores 2 --configfile "$CONFIG" --until blastn
        local blastn="test/results/test/test_R/BLASTN/test_R.contigs.blastn.csv"
        [[ -s "$blastn" ]] || { echo "[smoke] FAIL: missing $blastn"; exit 1; }
        echo "[smoke] OK: $blastn"
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
