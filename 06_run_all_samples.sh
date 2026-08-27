#!/bin/bash
# ============================================================
# 06_run_all_samples.sh
# Submits STARsolo + CellBender for all 28 samples (GSE204770).
# Can be called directly or from 05_run_pipeline.sh.
#
# Usage:
#   bash 06_run_all_samples.sh                        # run immediately
#   bash 06_run_all_samples.sh --dependency=afterok:123  # wait for job 123
#
# Each sample gets its own pair of jobs:
#   starsolo → cellbender (cellbender waits for starsolo to finish)
# ============================================================

set -euo pipefail
mkdir -p logs

# Optional dependency flag passed in from 05_run_pipeline.sh
INDEX_DEP="${1:-}"

SAMPLES=(
    SRR19391326
    SRR19391327
    SRR19391328
    SRR19391329
    SRR19391330
    SRR19391331
    SRR19391332
    SRR19391333
    SRR19391334
    SRR19391335
    SRR19391336
    SRR19391337
    SRR19394872
    SRR19394876
    SRR19394877
    SRR19394878
    SRR19394879
    SRR19394880
    SRR19394886
    SRR19394887
    SRR19394888
    SRR19394889
    SRR19394890
    SRR19394891
    SRR19394896
    SRR19394897
    SRR19394898
    SRR19394899
)

echo "Submitting jobs for ${#SAMPLES[@]} samples..."
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    # STARsolo — waits for index if dependency passed in
    JID_STAR=$(sbatch \
        $INDEX_DEP \
        --job-name="star_${SAMPLE}" \
        --export=ALL,SAMPLE_NAME="$SAMPLE" \
        03_starsolo.sh | awk '{print $NF}')

    # CellBender — waits for STARsolo to finish
    JID_CB=$(sbatch \
        --dependency=afterok:"$JID_STAR" \
        --job-name="cb_${SAMPLE}" \
        --export=ALL,SAMPLE_NAME="$SAMPLE" \
        04_cellbender.sh | awk '{print $NF}')

    echo "$SAMPLE  →  STARsolo: $JID_STAR  |  CellBender: $JID_CB"
done

echo ""
echo "All jobs submitted."
echo "Monitor with: squeue -u \$USER"
echo "Watch live:   watch -n 30 squeue -u \$USER"
