#!/bin/bash
#SBATCH --job-name=cellbender.SBATCH
#SBATCH --output=z04.%x.%j
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrb339@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=90:00:00
#SBATCH --mem=120G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# ============================================================
# 04_cellbender.sh
# CellBender remove-background for all 26 samples (GSE204770).
# Input: Cell Ranger raw h5 files from GEO.
# Submit with: sbatch 04_cellbender.sh
#
# If your cluster has no GPU partition, remove --partition=gpu,
# --gres=gpu:1, and --cuda below — expect much longer runtime.
# ============================================================

set -euo pipefail
mkdir -p logs

CELLBENDER_BIN=~/.conda/envs/cellbender_env/bin/cellbender

# All 26 samples from GSE204770_RAW.tar
SAMPLES=(
    "GSM6190636_A1"
    "GSM6190637_A2"
    "GSM6190638_A3"
    "GSM6190639_A4"
    "GSM6190640_A5"
    "GSM6190641_B7"
    "GSM6190642_B8"
    "GSM6190643_B9"
    "GSM6190829_1"
    "GSM6190830_2"
    "GSM6190831_3"
    "GSM6190832_4"
    "GSM6190833_5"
    "GSM6190834_6"
    "GSM6190835_7"
    "GSM6190836_8"
    "GSM6190837_9"
    "GSM6190838_10"
    "GSM6190839_11"
    "GSM6190840_12"
    "GSM6190841_13"
    "GSM6190842_14"
    "GSM6190843_15"
    "GSM6190844_16"
    "GSM6190845_17"
    "GSM6190846_18"
)

INPUT_DIR="$HOME/geo_raw"
OUT_BASE="$HOME/cellbender_out"

# Adjust EXPECTED_CELLS based on the paper's reported nuclei counts.
# Paper reported 58,079 nuclei across all samples, roughly ~2,000-3,000 per sample.
EXPECTED_CELLS=3000
TOTAL_DROPLETS=15000
FPR="0.01"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=== CellBender: $SAMPLE ==="

    INPUT_H5="$INPUT_DIR/${SAMPLE}_raw_feature_bc_matrix.h5"
    OUT_DIR="$OUT_BASE/$SAMPLE"
    OUT_FILE="$OUT_DIR/${SAMPLE}_cellbender.h5"

    if [[ ! -f "$INPUT_H5" ]]; then
        echo "ERROR: Input file not found at $INPUT_H5 — skipping"
        continue
    fi

    mkdir -p "$OUT_DIR"

    $CELLBENDER_BIN remove-background \
        --input "$INPUT_H5" \
        --output "$OUT_FILE" \
        --expected-cells "$EXPECTED_CELLS" \
        --total-droplets-included "$TOTAL_DROPLETS" \
        --fpr "$FPR" \
        --epochs 150 \
        --cuda \
        2>&1 | tee "$OUT_DIR/cellbender_run.log"

    echo "=== Done: $SAMPLE ==="
    echo "Cleaned h5: ${OUT_FILE%.h5}_filtered.h5"
done

echo ""
echo "=== All samples complete ==="
