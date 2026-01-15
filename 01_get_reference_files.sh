#!/bin/bash
set -euo pipefail

# Download mm10 genome FASTA (UCSC) and GENCODE M25 annotation (GTF).
# Run this from the directory where you want the files saved (e.g., 01.RawData).

MM10_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

echo "Downloading mm10 FASTA..."
wget -O mm10.fa.gz "${MM10_URL}"
gunzip -f mm10.fa.gz

echo "Downloading GENCODE M25 GTF..."
wget -O gencode.vM25.annotation.gtf.gz "${GTF_URL}"
gunzip -f gencode.vM25.annotation.gtf.gz

echo "Done. Created: mm10.fa and gencode.vM25.annotation.gtf"

