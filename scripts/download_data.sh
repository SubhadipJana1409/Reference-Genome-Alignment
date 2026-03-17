#!/usr/bin/env bash
# =============================================================================
# download_data.sh — Download reference genome and annotation for Day 2
#
# Downloads:
#   - GRCh38 chr22 FASTA  (Ensembl release 110)
#   - GENCODE v44 chr22 GTF annotation
#
# Usage:
#   bash scripts/download_data.sh
#   bash scripts/download_data.sh --synthetic   # generate synthetic chr22 instead
#
# #30DaysOfBioinformatics | SubhadipJana1409
# =============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="$PROJECT_DIR/data/reference"
WGS_DIR="$PROJECT_DIR/data/wgs"
RNA_DIR="$PROJECT_DIR/data/rnaseq"

mkdir -p "$REF_DIR" "$WGS_DIR" "$RNA_DIR"

# ── Parse arguments ──────────────────────────────────────────────────────────
SYNTHETIC=false
for arg in "$@"; do
    [[ "$arg" == "--synthetic" ]] && SYNTHETIC=true
done

echo "================================================"
echo " Day 2 — Reference Data Download"
echo "================================================"
echo " Output: $REF_DIR"
echo ""

# ── Reference genome ─────────────────────────────────────────────────────────
if [[ "$SYNTHETIC" == "true" ]]; then
    echo "[1/2] Generating synthetic chr22 (seed=42, 5 Mb)..."
    python3 scripts/simulate_reference.py
    echo "      → $REF_DIR/chr22.fa"
else
    echo "[1/2] Downloading GRCh38 chr22 from Ensembl release 110..."
    CHR22_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
    CHR22_GZ="$REF_DIR/chr22.fa.gz"

    if [[ -f "$REF_DIR/chr22.fa" ]]; then
        echo "      Already exists — skipping download."
    else
        wget -q --show-progress -O "$CHR22_GZ" "$CHR22_URL"
        echo "      Decompressing..."
        gunzip "$CHR22_GZ"
        # Rename header to match pipeline expectation
        sed -i '1s/.*/\>chr22/' "$REF_DIR/chr22.fa"
        echo "      → $REF_DIR/chr22.fa ($(du -sh "$REF_DIR/chr22.fa" | cut -f1))"
    fi
fi

# ── GENCODE v44 chr22 GTF ─────────────────────────────────────────────────────
echo ""
echo "[2/2] Downloading GENCODE v44 chr22 annotation..."

GTF_OUT="$REF_DIR/gencode_v44_chr22_gtf"
FULL_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
FULL_GTF_GZ="/tmp/gencode.v44.annotation.gtf.gz"

if [[ -f "$GTF_OUT" ]]; then
    echo "      Already exists — skipping download."
else
    # Download full GTF then extract chr22 (more reliable than chr22-specific URL)
    if [[ ! -f "$FULL_GTF_GZ" ]]; then
        echo "      Downloading full GENCODE v44 GTF (~500 MB)..."
        wget -q --show-progress -O "$FULL_GTF_GZ" "$FULL_GTF_URL"
    fi

    echo "      Extracting chr22 annotations..."
    zcat "$FULL_GTF_GZ" | awk '$1=="chr22" || /^#/' > "$GTF_OUT"
    echo "      → $GTF_OUT ($(wc -l < "$GTF_OUT") lines)"

    # Clean up full GTF to save space
    rm -f "$FULL_GTF_GZ"
fi

# ── Simulate reads with wgsim ─────────────────────────────────────────────────
echo ""
echo "[3/3] Simulating reads with wgsim..."

if command -v wgsim &>/dev/null; then
    REF="$REF_DIR/chr22.fa"

    # WGS samples
    echo "  WGS (NA12878) — 3 samples × 100k pairs × 150bp..."
    declare -A WGS_PARAMS=(
        ["S1"]="0.01 400 42"
        ["S2"]="0.015 400 100"
        ["S3"]="0.01 350 200"
    )
    for s in S1 S2 S3; do
        read ERR INSERT SEED <<< "${WGS_PARAMS[$s]}"
        R1="$WGS_DIR/NA12878_${s}_R1.fastq"
        R2="$WGS_DIR/NA12878_${s}_R2.fastq"
        if [[ ! -f "${R1}.gz" ]]; then
            wgsim -N 100000 -1 150 -2 150 -e "$ERR" -d "$INSERT" -s 50 -S "$SEED" \
                  "$REF" "$R1" "$R2" > /dev/null 2>&1
            gzip "$R1" "$R2"
            echo "    NA12878_${s}: done"
        else
            echo "    NA12878_${s}: already exists"
        fi
    done

    # RNA-seq samples
    echo "  RNA-seq (PBMC) — 3 samples × 100k pairs × 150bp..."
    declare -A RNA_PARAMS=(
        ["S1"]="0.005 300 42"
        ["S2"]="0.005 300 100"
        ["S3"]="0.008 280 200"
    )
    for s in S1 S2 S3; do
        read ERR INSERT SEED <<< "${RNA_PARAMS[$s]}"
        R1="$RNA_DIR/PBMC_${s}_R1.fastq"
        R2="$RNA_DIR/PBMC_${s}_R2.fastq"
        if [[ ! -f "${R1}.gz" ]]; then
            wgsim -N 100000 -1 150 -2 150 -e "$ERR" -d "$INSERT" -s 50 -S "$SEED" \
                  "$REF" "$R1" "$R2" > /dev/null 2>&1
            gzip "$R1" "$R2"
            echo "    PBMC_${s}: done"
        else
            echo "    PBMC_${s}: already exists"
        fi
    done
else
    echo "  wgsim not found — install samtools (includes wgsim) and re-run."
    echo "  Or provide your own FASTQ files in data/wgs/ and data/rnaseq/"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "================================================"
echo " Done! Files ready:"
echo "================================================"
find "$REF_DIR" "$WGS_DIR" "$RNA_DIR" -type f | sort | while read f; do
    printf "  %-60s %s\n" "${f#$PROJECT_DIR/}" "$(du -sh "$f" | cut -f1)"
done

echo ""
echo "Next step:"
echo "  python scripts/alignment_benchmark.py --index-only   # build indexes"
echo "  python scripts/alignment_benchmark.py                # full pipeline"
