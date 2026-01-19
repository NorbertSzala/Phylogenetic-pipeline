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

SPECIES="$1"
OUTPUT="$2"

if [[ -z "${SPECIES}" || -z "${OUTPUT}" ]]; then
    echo "Usage: select_genomes.sh <species_name> <output_json>" >&2
    exit 1
fi

echo "species=${SPECIES}"

TMP=$(mktemp)

# Run NCBI datasets
if ! datasets summary genome taxon "${SPECIES}" > "${TMP}"; then
    echo "[select_genomes] datasets failed for ${SPECIES}" >&2
    rm -f "${TMP}"
    exit 2
fi

# Check non-empty JSON
if [[ ! -s "${TMP}" ]]; then
    echo "[select_genomes] empty result for ${SPECIES}" >&2
    rm -f "${TMP}"
    exit 3
fi

mv "${TMP}" "${OUTPUT}"

echo "[select_genomes] summary written to ${OUTPUT}"
