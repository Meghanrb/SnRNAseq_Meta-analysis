#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=standard

# ============================================================
# 02_download_fastq.sh
# Downloads FASTQs from SRA/GEO using prefetch + fasterq-dump.
# Each SRR accession in SRR_LIST gets its own subdirectory.
#
# EDIT the variables in the USER CONFIG block before submitting.
# Submit once per sample (or loop over SRR_LIST).
# ============================================================

set -euo pipefail
mkdir -p logs

# ---- USER CONFIG -------------------------------------------
# Space-separated list of SRR accessions for ONE 10x sample.
# A single 10x run can have multiple SRR accessions (lanes).
# R1 = barcode+UMI (28 bp), R2 = cDNA read.
SRR_LIST="SRR0000001 SRR0000002"    # <-- REPLACE with your accessions
SAMPLE_NAME="sample1"               # used to name output directories
FASTQ_OUT="$HOME/fastqs/$SAMPLE_NAME"
THREADS=${SLURM_CPUS_PER_TASK:-8}
TMP_DIR="$HOME/sra_cache"
# ---- END USER CONFIG ---------------------------------------

CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate star_env

mkdir -p "$FASTQ_OUT" "$TMP_DIR"

# Point sra-tools at a writable cache
vdb-config --set /repository/user/main/public/root="$TMP_DIR" 2>/dev/null || true

for SRR in $SRR_LIST; do
    echo "=== Downloading $SRR ==="
    prefetch --output-directory "$TMP_DIR" "$SRR"
    fasterq-dump \
        --outdir "$FASTQ_OUT" \
        --threads "$THREADS" \
        --split-files \
        --temp "$TMP_DIR" \
        "$TMP_DIR/$SRR/$SRR.sra"
    # Compress to save space
    pigz -p "$THREADS" "$FASTQ_OUT/${SRR}"_*.fastq
    echo "$SRR done."
done

echo "FASTQs in: $FASTQ_OUT"
conda deactivate
