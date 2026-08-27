#!/bin/bash
#SBATCH --job-name=install_envs
#SBATCH --output=logs/install_envs_%j.out
#SBATCH --error=logs/install_envs_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=standard

# ============================================================
# 00_install_envs.sh
# Run ONCE to install conda environments for STARsolo and CellBender.
# Edit CONDA_BASE if your conda/mamba lives elsewhere.
# ============================================================

set -euo pipefail
mkdir -p logs

# ---- locate conda ------------------------------------------
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Prefer mamba for speed; fall back to conda
if command -v mamba &>/dev/null; then
    PKG=mamba
else
    PKG=conda
fi

echo "Using: $PKG  |  conda base: $CONDA_BASE"

# ---- STARsolo environment ----------------------------------
echo "=== Creating star_env ==="
$PKG create -y -n star_env -c conda-forge -c bioconda \
    star=2.7.11a \
    samtools=1.20 \
    sra-tools=3.1.1 \
    pigz \
    python=3.11

# ---- CellBender environment (CPU-only fallback) ------------
# If your cluster has GPUs, the GPU-aware PyTorch below is preferred.
# For CPU-only clusters, replace the torch line with:
#   pytorch cpuonly -c pytorch
echo "=== Creating cellbender_env ==="
$PKG create -y -n cellbender_env -c conda-forge -c pytorch \
    python=3.10 \
    pytorch pytorch-cuda=11.8 -c pytorch -c nvidia

conda activate cellbender_env
pip install cellbender
conda deactivate

echo "=== All environments installed ==="
echo "STARsolo env : star_env"
echo "CellBender env: cellbender_env"
