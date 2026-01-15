#!/bin/bash
#SBATCH --job-name=02_STAR_Index_mm10
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load star-2.7.10b-gcc-12.1.0

# Adjust these paths if you run from somewhere else.
WORKDIR="/scratch/${USER}/rnaseq_project/01.RawData"
GENOME_FA="${WORKDIR}/mm10.fa"
GTF="${WORKDIR}/gencode.vM25.annotation.gtf"
STAR_INDEX="${WORKDIR}/STAR_index_mm10"

mkdir -p "${STAR_INDEX}"

# For 2x100bp reads, sjdbOverhang should be readLength-1 = 99, check your Novogene submission specifics
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "${STAR_INDEX}" --genomeFastaFiles "${GENOME_FA}" --sjdbGTFfile "${GTF}" --sjdbOverhang 99

echo "STAR index generated at: ${STAR_INDEX}"

