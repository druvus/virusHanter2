#!/usr/bin/env bash
#
# Smoke runner for virusHanter2.
#
#   ./test/run_smoke.sh           # lint + dry-run only (no reference data)
#   ./test/run_smoke.sh --full    # full pipeline against mini databases
#
# The dry-run path exercises DAG construction, rule imports, and the
# generate_report shell-out to `reporthanter`; it does not require any
# reference databases. Use --full only after populating test/mini_db/ as
# described in test/README.md.

set -euo pipefail

cd "$(dirname "$0")/.."

CONFIG="test/config.test.yaml"

echo "[smoke] snakemake --lint --configfile $CONFIG"
snakemake --lint --configfile "$CONFIG"

echo "[smoke] snakemake -n --sdm conda --configfile $CONFIG"
snakemake -n --sdm conda --configfile "$CONFIG"

if [[ "${1:-}" == "--full" ]]; then
    echo "[smoke] snakemake --sdm conda --cores 2 --configfile $CONFIG"
    snakemake --sdm conda --cores 2 --configfile "$CONFIG"

    sample_dir="test/results/test"
    html="$sample_dir/test/REPORT/test.html"
    run_csv="$sample_dir/run_information_test.csv"

    if [[ ! -s "$html" ]]; then
        echo "[smoke] FAIL: expected $html"
        exit 1
    fi
    if [[ ! -s "$run_csv" ]]; then
        echo "[smoke] FAIL: expected $run_csv"
        exit 1
    fi
    echo "[smoke] OK: $html, $run_csv"
fi
