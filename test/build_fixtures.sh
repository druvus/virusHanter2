#!/usr/bin/env bash
#
# Thin wrapper around the fixture-build Snakemake workflow.
# All actual fixture rules live in test/build_fixtures.smk; this script
# just forwards to Snakemake so the per-tool conda envs (envs/*.yaml) are
# materialised automatically.
#
# Required on PATH (driver env): snakemake>=9, conda (or mamba).
# Every bioinformatics tool — bwa, kraken2-build, kaiju-mkbwt, kaiju-mkfmi,
# makeblastdb — is supplied by the per-rule conda env from envs/.

set -euo pipefail

cd "$(dirname "$0")/.."

CORES="${SMOKE_CORES:-2}"

# osx-arm64 fallback (see run_smoke.sh for context).
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    export CONDA_SUBDIR="${CONDA_SUBDIR:-osx-64}"
fi

snakemake \
    --snakefile test/build_fixtures.smk \
    --sdm conda \
    --cores "$CORES" \
    --directory .

echo "[fixtures] done. Tree:"
find test/mini_db -maxdepth 2 -print
