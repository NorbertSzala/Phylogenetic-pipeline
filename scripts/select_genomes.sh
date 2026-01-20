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

# Cache
if [[ -s "${OUTPUT}" ]]; then
    echo "[select_genomes] cached: ${OUTPUT}"
    exit 0
fi

if ! datasets summary genome taxon "${SPECIES}" \
    --assembly-level complete \
    --reference \
    --limit 1 \
    > "${OUTPUT}"; then
    echo "[select_genomes] datasets failed for ${SPECIES}, creating placeholder" >&2
    : > "${OUTPUT}"        # <<< kluczowe dla Nextflow
    exit 0
fi

if [[ ! -s "${OUTPUT}" ]]; then
    echo "[select_genomes] empty result for ${SPECIES}, keeping placeholder" >&2
    : > "${OUTPUT}"
    exit 0
fi

echo "[select_genomes] summary written to ${OUTPUT}"



