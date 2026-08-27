#!/bin/bash
# ============================================================
# 05_run_pipeline.sh
# Submits all pipeline steps in dependency order for GSE204770.
# Run this from the login node: bash 05_run_pipeline.sh
#
# Steps:
#   00 → install conda environments (run once; skip with --skip-install)
#   01 → build STAR index for mm10  (run once; skip with --skip-index)
#   03 → STARsolo alignment for all 28 samples
#   04 → CellBender for all 28 samples (depends on 03)
#
# Prerequisites before running:
#   - FASTQs transferred to /home/mrb339/downloaded_data_hou/<SRR>/
#   - Barcode whitelist at $HOME/references/whitelists/737K-august-2016.txt
# ============================================================

set -euo pipefail
mkdir -p logs

SKIP_INSTALL=false
SKIP_INDEX=false

for arg in "$@"; do
    [[ "$arg" == "--skip-install" ]] && SKIP_INSTALL=true
    [[ "$arg" == "--skip-index"   ]] && SKIP_INDEX=true
done

# Step 0: install environments
if $SKIP_INSTALL; then
    echo "Skipping 00_install_envs.sh (--skip-install)"
    JID0=""
else
    JID0=$(sbatch 00_install_envs.sh | awk '{print $NF}')
    echo "Submitted 00_install_envs  → job $JID0"
fi

# Step 1: build STAR index (depends on install)
if $SKIP_INDEX; then
    echo "Skipping 01_build_star_index.sh (--skip-index)"
    JID1=""
else
    DEP="${JID0:+--dependency=afterok:$JID0}"
    JID1=$(sbatch $DEP 01_build_star_index.sh | awk '{print $NF}')
    echo "Submitted 01_build_star_index → job $JID1"
fi

# Step 2: STARsolo for all 28 samples (depends on index)
DEP="${JID1:+--dependency=afterok:$JID1}"
JID2=$(sbatch $DEP 03_starsolo.sh | awk '{print $NF}')
echo "Submitted 03_starsolo         → job $JID2"

# Step 3: CellBender for all 28 samples (depends on STARsolo)
JID3=$(sbatch --dependency=afterok:"$JID2" 04_cellbender.sh | awk '{print $NF}')
echo "Submitted 04_cellbender       → job $JID3"

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Watch live: watch -n 30 squeue -u \$USER"
