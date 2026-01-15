#!/bin/bash
set -euo pipefail

# Usage:
# bash 04_install_subread_featurecounts.sh /scratch/$USER/tools

INSTALL_BASE="${1:-/scratch/${USER}/tools}"
mkdir -p "${INSTALL_BASE}" # makes new directory if install base name directory is not found
cd "${INSTALL_BASE}"

URL="https://downloads.sourceforge.net/project/subread/subread-2.1.0/subread-2.1.0-Linux-x86_64.tar.gz"
ARCHIVE="subread-2.1.0-Linux-x86_64.tar.gz"
DIR="subread-2.1.0-Linux-x86_64"

echo "Downloading Subread..."
wget -O "${ARCHIVE}" "${URL}"

echo "Extracting..."
tar -xzf "${ARCHIVE}" # must use -xzf flags coming from *.tar.gz

echo "Installed at: ${INSTALL_BASE}/${DIR}"
echo "featureCounts binary: ${INSTALL_BASE}/${DIR}/bin/featureCounts"

