#!/bin/bash
#SBATCH --job-name=05_featureCounts
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# ---- EDIT THESE PATHS ----
PROJECT_DIR="/scratch/${USER}/rnaseq_project"
GTF="${PROJECT_DIR}/01.RawData/gencode.vM25.annotation.gtf"
SUBREAD_DIR="/scratch/${USER}/tools/subread-2.1.0-Linux-x86_64"
# --------------------------

export PATH="${SUBREAD_DIR}/bin:${PATH}"

# Example: BAMs produced by Step 3 are under alignments/<sample>/<sample>_Aligned.sortedByCoord.out.bam

BAMS=(
  "${PROJECT_DIR}/alignments/sample1/sample1_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/A2/A2_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/B1/B1_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/B2/B2_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/B3/B3_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/Mix_A3_4674_7281/Mix_A3_4674_7281_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/Mix_A4_4675_7282/Mix_A4_4675_7282_Aligned.sortedByCoord.out.bam"
  "${PROJECT_DIR}/alignments/Mix_B4_4675_7286/Mix_B4_4675_7286_Aligned.sortedByCoord.out.bam"
)

cd "${PROJECT_DIR}"

# -T uses threads (match cpus-per-task)
# -p indicates paired-end
# -B requires both ends mapped
# -C checks chimeric fragments
featureCounts -T 8 -p -B -C -a "${GTF}" -o all_counts.txt "${BAMS[@]}"

echo "featureCounts complete: ${PROJECT_DIR}/all_counts.txt"

