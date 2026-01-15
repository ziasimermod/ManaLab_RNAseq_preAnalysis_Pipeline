# RNA-seq Pre-analysis Pipeline (QC → trim → align → counts)

### This repository contains a simple RNA-seq pre-analysis workflow intended for use on an HPC cluster running SLURM. It performs:

#### Download reference genome + annotation (mm10 + GENCODE M25)

#### Build STAR genome index

#### Run per-sample FastQC → Trimmomatic → STAR alignment (as a SLURM array)

#### Install Subread (featureCounts)

#### Generate a gene-level count matrix with featureCounts

## Requirements
#### Cluster / Scheduler

#### SLURM

#### Available environment modules

## Load modules using:

module load <MODULE_NAME>


## Modules assumed available:

fastqc-0.12.1-gcc-11.2.0

trimmomatic-0.39-gcc-12.1.0

star-2.7.10b-gcc-12.1.0

samtools-1.21-gcc-12.1.0 (required for BAM indexing)

If the cluster uses slightly different module names, adjust the module load lines in the scripts. Use 'module avail' to see available modules

## Input data expectations

This pipeline assumes paired-end FASTQ files exist and follow one of these patterns:

#### Option A (per-sample subfolders):

/path/to/RawData/<SAMPLE_ID>/<SAMPLE_ID>_1.fastq.gz
/path/to/RawData/<SAMPLE_ID>/<SAMPLE_ID>_2.fastq.gz


#### Option B (flat folder):

/path/to/RawData/<SAMPLE_ID>_1.fastq.gz
/path/to/RawData/<SAMPLE_ID>_2.fastq.gz


The provided array script is written for Option A. If you use Option B, modify how SAMPLE_R1 and SAMPLE_R2 are defined.

## Reference files

### This pipeline uses:

#### mm10 genome FASTA from UCSC

#### GENCODE mouse annotation GTF release M25

#### Downloaded by: scripts/01_get_reference_files.sh

# Step-by-step usage
## Step 0: choose a working directory

#### Example:

mkdir -p /scratch/$USER/rnaseq_project
cd /scratch/$USER/rnaseq_project
mkdir -p 01.RawData 02.Results logs scripts


#### Copy this repo’s scripts/ folder into your project or symlink it.

## Step 1: Download reference genome + GTF

From your project directory:

cd 01.RawData
bash ../scripts/01_get_reference_files.sh


### This creates:

#### mm10.fa

#### gencode.vM25.annotation.gtf

## Step 2: Build STAR index

Edit scripts/02_build_star_index.sh to confirm paths and threads, then run:

cd 01.RawData
sbatch ../scripts/02_build_star_index.sh


### Monitor jobs with:

sacct
squeue -u $USER

## Step 3: Run QC → trim → align (SLURM array)

Edit scripts/03_rnaseq_qc_trim_align_array.sh:

set SAMPLES_DIR

set STAR_INDEX

update the samples=(...) list

set --array=0-(N-1) for your number of samples

Then submit:

sbatch ../scripts/03_rnaseq_qc_trim_align_array.sh


Outputs include:

fastqc_results/<sample>/

trimmed FASTQs

STAR-aligned BAMs + indices

STAR gene counts (optional output from STAR)

## Step 4: Install Subread (featureCounts)

Run once (interactive login node is fine if allowed; otherwise use a small SLURM job):

bash scripts/04_install_subread_featurecounts.sh /scratch/$USER/tools


This installs to:

/scratch/$USER/tools/subread-2.1.0-Linux-x86_64/

## Step 5: Generate gene-level counts with featureCounts

Edit scripts/05_run_featurecounts.sh:

set SUBREAD_DIR

set GTF

list your BAM files (or use a glob if your naming is consistent)

Submit:

sbatch scripts/05_run_featurecounts.sh


## Output:

#### all_counts.txt (featureCounts output table)

Notes / common debugs

Do not modify .bashrc inside batch jobs. Use export PATH=...:$PATH in the job script (already done below).

Avoid multi-line \ in scripts if your environment occasionally mis-parses them. All commands here are written as single-line invocations.

Ensure your FASTQ naming matches ***exactly*** what the scripts expect.

STAR indexing can take time; be sure you allocate enough memory/time for Step 2 on your cluster.
