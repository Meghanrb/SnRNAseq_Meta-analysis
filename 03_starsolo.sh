#!/bin/bash
#SBATCH --job-name=starsolo.SBATCH
#SBATCH --output=z03.%x.%j
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrb339@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=90:00:00
#SBATCH --mem=120G

# ============================================================
# 03_starsolo.sh
# STARsolo alignment for all 17 samples (GSE204770).
# 10x Chromium 5' v1 snRNA-seq, mouse mm10.
# Multi-lane samples: both SRRs are fed into STAR together.
# Submit with: sbatch 03_starsolo.sh
# ============================================================

set -euo pipefail
mkdir -p logs

STAR_BIN=~/.conda/envs/star_env/bin/STAR
SAMTOOLS_BIN=~/.conda/envs/star_env/bin/samtools

# Format: "GSM_ID:SRR1,SRR2" or "GSM_ID:SRR1" for single-lane samples
SAMPLES=(
    "GSM6190643:SRR19391326,SRR19391327"
    "GSM6190642:SRR19391328,SRR19391329"
    "GSM6190641:SRR19391330,SRR19391331"
    "GSM6190637:SRR19391335,SRR19391336"
    "GSM6190840:SRR19394877,SRR19394878"
    "GSM6190839:SRR19394879,SRR19394880"
    "GSM6190835:SRR19394886,SRR19394887"
    "GSM6190834:SRR19394888,SRR19394889"
    "GSM6190833:SRR19394890,SRR19394891"
    "GSM6190830:SRR19394896,SRR19394897"
    "GSM6190829:SRR19394898,SRR19394899"
    "GSM6190640:SRR19391332"
    "GSM6190639:SRR19391333"
    "GSM6190638:SRR19391334"
    "GSM6190636:SRR19391337"
    "GSM6190845:SRR19394872"
    "GSM6190841:SRR19394876"
)

FASTQ_BASE="/home/mrb339/downloaded_data_hou"
GENOME_DIR="$HOME/references/GRCm38_star_index"
WHITELIST="$HOME/references/whitelists/737K-august-2016.txt"
THREADS=${SLURM_CPUS_PER_TASK:-1}
CB_LEN=16
UMI_LEN=10

echo "STAR binary: $STAR_BIN"
echo "STAR version: $($STAR_BIN --version)"

if [[ ! -f "$WHITELIST" ]]; then
    echo "ERROR: Whitelist not found at $WHITELIST"
    exit 1
fi

for ENTRY in "${SAMPLES[@]}"; do
    GSM="${ENTRY%%:*}"
    SRRS="${ENTRY##*:}"

    echo ""
    echo "=== Processing $GSM (SRRs: $SRRS) ==="

    OUT_DIR="$HOME/starsolo_out/$GSM"
    mkdir -p "$OUT_DIR"

    # Collect R1 and R2 files from all SRR directories for this sample
    R1_FILES=""
    R2_FILES=""
    for SRR in ${SRRS//,/ }; do
        FASTQ_DIR="$FASTQ_BASE/$SRR"
        R1=$(ls "$FASTQ_DIR"/*_R1_001.fastq.gz 2>/dev/null | tr '\n' ',' | sed 's/,$//')
        R2=$(ls "$FASTQ_DIR"/*_R2_001.fastq.gz 2>/dev/null | tr '\n' ',' | sed 's/,$//')
        if [[ -z "$R1" || -z "$R2" ]]; then
            echo "ERROR: Could not find FASTQs in $FASTQ_DIR — skipping $GSM"
            continue 2
        fi
        R1_FILES="${R1_FILES:+$R1_FILES,}$R1"
        R2_FILES="${R2_FILES:+$R2_FILES,}$R2"
    done

    echo "R1: $R1_FILES"
    echo "R2: $R2_FILES"

    $STAR_BIN \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBstart 1 --soloCBlen "$CB_LEN" \
        --soloUMIstart 17 --soloUMIlen "$UMI_LEN" \
        --soloBarcodeReadLength 0 \
        \
        --readFilesIn "$R2_FILES" "$R1_FILES" \
        --readFilesCommand zcat \
        \
        --genomeDir "$GENOME_DIR" \
        --runThreadN "$THREADS" \
        \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        \
        --soloFeatures Gene GeneFull \
        --outFileNamePrefix "$OUT_DIR/" \
        \
        --soloCellFilter EmptyDrops_CR \
        \
        --outFilterScoreMin 30 \
        \
        --runDirPerm All_RWX \
        2>&1 | tee "$OUT_DIR/starsolo_run.log"

    $SAMTOOLS_BIN index -@ "$THREADS" "$OUT_DIR/Aligned.sortedByCoord.out.bam"

    echo "=== Done: $GSM ==="
    echo "Raw matrix:      $OUT_DIR/Solo.out/GeneFull/raw/"
    echo "Filtered matrix: $OUT_DIR/Solo.out/GeneFull/filtered/"
done

echo ""
echo "=== All samples complete ==="
