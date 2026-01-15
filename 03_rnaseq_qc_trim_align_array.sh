#!/bin/bash
#SBATCH --job-name=03_RNAseq_QC_Trim_Align
#SBATCH --array=0-7 # Adjust this based on the number of samples you are running, be sure to update your 'samples' too
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

module load fastqc-0.12.1-gcc-11.2.0
module load trimmomatic-0.39-gcc-12.1.0
module load star-2.7.10b-gcc-12.1.0
module load samtools-1.21-gcc-12.1.0

# ---- EDIT THESE PATHS ----
PROJECT_DIR="/scratch/${USER}/rnaseq_project"
SAMPLES_DIR="${PROJECT_DIR}/01.RawData"
STAR_INDEX="${SAMPLES_DIR}/STAR_index_mm10"
# --------------------------

# Sample IDs must match folder names if using /SAMPLES_DIR/<sample>/<sample>_1.fastq.gz, these need to be specified per pipeline run
samples=("A1" "A2" "B1" "B2" "B3" "Mix_A3_4674_7281" "Mix_A4_4675_7282" "Mix_B4_4675_7286")
sample="${samples[$SLURM_ARRAY_TASK_ID]}"

echo "Processing sample: ${sample}"

SAMPLE_R1="${SAMPLES_DIR}/${sample}/${sample}_1.fastq.gz"
SAMPLE_R2="${SAMPLES_DIR}/${sample}/${sample}_2.fastq.gz"

mkdir -p "${PROJECT_DIR}/fastqc_results/${sample}"
mkdir -p "${PROJECT_DIR}/alignments/${sample}"
mkdir -p "${PROJECT_DIR}/trimmed/${sample}"

# Step 1: FastQC
fastqc "${SAMPLE_R1}" "${SAMPLE_R2}" -o "${PROJECT_DIR}/fastqc_results/${sample}"

# Step 2: Trimmomatic (paired-end)
# NOTE: TruSeq3-PE.fa must be discoverable by Trimmomatic in your environment.
# If your cluster doesn’t provide it automatically, set an explicit path to the adapter FASTA.
TRIM_R1_PAIRED="${PROJECT_DIR}/trimmed/${sample}/${sample}_1_trimmed.fastq.gz"
TRIM_R1_UNPAIRED="${PROJECT_DIR}/trimmed/${sample}/${sample}_1_unpaired.fastq.gz"
TRIM_R2_PAIRED="${PROJECT_DIR}/trimmed/${sample}/${sample}_2_trimmed.fastq.gz"
TRIM_R2_UNPAIRED="${PROJECT_DIR}/trimmed/${sample}/${sample}_2_unpaired.fastq.gz"

trimmomatic PE -threads 4 "${SAMPLE_R1}" "${SAMPLE_R2}" "${TRIM_R1_PAIRED}" "${TRIM_R1_UNPAIRED}" "${TRIM_R2_PAIRED}" "${TRIM_R2_UNPAIRED}" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: STAR alignment
OUT_PREFIX="${PROJECT_DIR}/alignments/${sample}/${sample}_"
STAR --genomeDir "${STAR_INDEX}" --runThreadN 4 --readFilesIn "${TRIM_R1_PAIRED}" "${TRIM_R2_PAIRED}" --readFilesCommand zcat --outFileNamePrefix "${OUT_PREFIX}" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

# Step 4: Index BAM
BAM="${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
samtools index "${BAM}"

echo "Completed sample: ${sample}"

