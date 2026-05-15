#!/usr/bin/env bash
#
# Build the mini reference databases the smoke runner consumes.
#
# Tools required on PATH (a single conda env is fine):
#   python pandas pyarrow bwa samtools kraken2 kaiju makeblastdb
#
# CheckV is intentionally not built here. CheckV's database is large and
# not amenable to synthesis from a few kb of reference. The smoke runner
# degrades to BLASTN-only mode when no CheckV DB is present.

set -euo pipefail

cd "$(dirname "$0")/.."

MINI=test/mini_db
SCRIPTS=test/scripts

echo "[fixtures] cleaning previous artifacts"
rm -rf "$MINI"
mkdir -p "$MINI"

echo "[fixtures] synthesising FASTQ pair"
python "$SCRIPTS/synthesize_fastq.py" \
    --out-r1 test/test_R1.fastq \
    --out-r2 test/test_R2.fastq

echo "[fixtures] writing reference FASTA files"
python "$SCRIPTS/write_references.py" \
    --host-fasta  "$MINI/host.fasta" \
    --virus-fasta "$MINI/virus.fasta" \
    --virus-protein-fasta "$MINI/virus_aa.fasta"

# ------------------------------------------------------------------ BWA host
echo "[fixtures] building BWA host index"
mkdir -p "$MINI/human"
bwa index -p "$MINI/human/human" "$MINI/host.fasta" >/dev/null 2>&1

# ------------------------------------------------------------------ Kraken2
echo "[fixtures] building Kraken2 mini DB"
mkdir -p "$MINI/kraken/taxonomy" "$MINI/kraken/library"
# Minimal taxonomy: root -> Viruses -> synthetic_virus
cat > "$MINI/kraken/taxonomy/nodes.dmp" <<'EOF'
1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
10239	|	1	|	superkingdom	|	|	9	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
100001	|	10239	|	species	|	|	9	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|
EOF
cat > "$MINI/kraken/taxonomy/names.dmp" <<'EOF'
1	|	root	|		|	scientific name	|
10239	|	Viruses	|		|	scientific name	|
100001	|	synthetic virus	|		|	scientific name	|
EOF

# Tag the FASTA header for kraken2-build accelerated-from-library mode.
sed 's/^>synthetic_virus/>synthetic_virus|kraken:taxid|100001/' \
    "$MINI/virus.fasta" > "$MINI/kraken/library/virus.fna"

kraken2-build --add-to-library "$MINI/kraken/library/virus.fna" \
    --db "$MINI/kraken" --no-masking >/dev/null 2>&1
kraken2-build --build --db "$MINI/kraken" >/dev/null 2>&1

# ------------------------------------------------------------------ Kaiju
echo "[fixtures] building Kaiju mini DB"
mkdir -p "$MINI/kaiju"
# Kaiju needs an NCBI-style taxdump pair to look up taxids; reuse Kraken's.
cp "$MINI/kraken/taxonomy/nodes.dmp" "$MINI/kaiju/nodes.dmp"
cp "$MINI/kraken/taxonomy/names.dmp" "$MINI/kaiju/names.dmp"
# Kaiju expects protein headers with the taxid right after the accession.
sed 's/^>synthetic_virus_aa/>P00001_100001/' \
    "$MINI/virus_aa.fasta" > "$MINI/kaiju/virus_aa.faa"

(
    cd "$MINI/kaiju"
    kaiju-mkbwt -n 2 -o kaiju_db virus_aa.faa >/dev/null 2>&1
    kaiju-mkfmi kaiju_db >/dev/null 2>&1
    rm -f kaiju_db.bwt kaiju_db.sa
)

# ------------------------------------------------------------------ BLAST nt
echo "[fixtures] building BLAST nt mini DB"
mkdir -p "$MINI/blast"
makeblastdb -in "$MINI/virus.fasta" -dbtype nucl \
    -out "$MINI/blast/viral" -title "viral_mini" >/dev/null 2>&1

# ------------------------------------------------------------------ Parquet
echo "[fixtures] writing virus.parquet"
python "$SCRIPTS/make_virus_parquet.py" --out "$MINI/virus.parquet"

# ------------------------------------------------------------------ CheckV stub
# Empty directory. The smoke runner skips CheckV when this is the only thing
# present; if you populate it with a real CheckV DB the smoke runs the full
# pipeline including the HTML report.
mkdir -p "$MINI/checkv"
: > "$MINI/checkv/.stub"

echo "[fixtures] done. Tree:"
find "$MINI" -maxdepth 2 -print
