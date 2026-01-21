#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# select_genomes.sh
#
# Download genome assembly summary from NCBI for ONE species.
#
# Designed for Nextflow:
# - no manual logging
# - stdout: progress information
# - stderr: errors
# - exit code != 0 signals failure
#
# Usage:
#   select_genomes.sh "<species_name>" <output_json>
# -----------------------------------------------------------------------------

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
