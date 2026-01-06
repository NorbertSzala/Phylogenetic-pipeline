#!/usr/bin/env bash

# Script to select genome sequences from NCBI given a list of gene species.

# Returns the parameters for downloading genomes, like assembly, L50, N50 etc.
# 1st part of pipeline

#Input: list of species names in ../data/short_taxonomy.csv
# Output: ../data/proteomes/selection/files with selected assemblies for select the best assemblies in select_best_assemblies.py

# one specie -> one summary.json

# Next step: select_best_assemblies.py

echo "#1. Selecting genomes from NCBI based on species names..."

# Safety settings - stop script on errors, undefined variables, or failed pipes
set -euo pipefail

INPUT=$1
OUTPUT=$2

if [[ -z "${INPUT}" || -z "${OUTPUT}" ]]; then
    echo "Usage: select_genomes.sh <species_name> <output_json>"
    exit 1
fi

# CREATE temp file
TMP=$(mktemp)

# Set logging
LOG="selecting_assemblies.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"
}

err() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" | tee -a "${LOG}" >&2
}


log "Starting genome selection"
log "Species: ${INPUT}"
log "Output: ${OUTPUT}"


if datasets summary genome taxon "$INPUT" > "${TMP}"; then
    # save only when JSON is not empty
    if [[ -s "${TMP}" ]]; then
        mv "${TMP}" "${OUTPUT}"
        log "Summary saved to ${OUTPUT}"

    else
        err "Empty result returned for species: ${INPUT}"

        rm -f "${TMP}"
        exit 2
    fi
else
    err "datasets command failed for species: ${INPUT}"
    rm -f "${TMP}"
    exit 3
fi
