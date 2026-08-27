#!/bin/bash
#SBATCH --job-name=star_index.SBATCH
#SBATCH --output=z01.%x
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mrb339@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=90:00:00
#SBATCH --mem=120G

# ============================================================
# 01_build_star_index.sh
# Builds a STAR genome index for mouse mm10 (GRCm38).
# Run once. Takes ~2-4 hrs with 1 CPU.
#
# NOTE: The original paper added a custom Clec7a cDNA to their
# reference. If you need results directly comparable to theirs,
# you must add that sequence to FASTA/GTF before running this.
# Contact the authors for the custom sequence.
# ============================================================

set -euo pipefail
mkdir -p logs

GENOME_DIR="$HOME/references/GRCm38_star_index"
FASTA_URL="https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
THREADS=${SLURM_CPUS_PER_TASK:-1}

CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate star_env

WORK_DIR="$HOME/references/downloads"
mkdir -p "$WORK_DIR" "$GENOME_DIR"

echo "Downloading mm10 (GRCm38 release 102) reference files..."
wget -nc -P "$WORK_DIR" "$FASTA_URL"
wget -nc -P "$WORK_DIR" "$GTF_URL"

FASTA_GZ="$WORK_DIR/$(basename $FASTA_URL)"
GTF_GZ="$WORK_DIR/$(basename $GTF_URL)"

echo "Decompressing..."
gunzip -kf "$FASTA_GZ"
gunzip -kf "$GTF_GZ"

FASTA="${FASTA_GZ%.gz}"
GTF="${GTF_GZ%.gz}"

echo "Building STAR index (1 CPU, expect 2-4 hrs)..."
STAR \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 149

echo "Index written to: $GENOME_DIR"
conda deactivate
