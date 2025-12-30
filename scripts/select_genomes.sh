#!/usr/bin/env bash

# Script to select genome sequences from NCBI given a list of gene species.

# Returns the parameters for downloading genomes, like assembly, L50, N50 etc.
# 1st part of pipeline

#Input: list of species names in ../data/short_taxonomy.csv
# Output: ../data/proteomes/selection with selected assemblies for select the best assemblies in select_best_assemblies.py

# Next step: select_best_assemblies.py

echo "#1. Selecting genomes from NCBI based on species names..."

# Safety settings - stop script on errors, undefined variables, or failed pipes
set -euo pipefail

INPUT=$1
OUTPUT=$2

# Set output and input directories/files
mkdir -p "$OUTPUT" ./logs

# Set logging
LOG="selecting_assemblies.log"
: > "$LOG" # clear log file

# read species names line by line from the input file
while IFS= read -r SPECIES; do
    [[ -z "$SPECIES" ]] && continue

    SAFE=$(echo "$SPECIES" | tr ' /' '__')
    TMP=$(mktemp)

    if datasets summary genome taxon "$SPECIES" > "$TMP"; then
        mv "$TMP" "${OUTPUT}/${SAFE}_summary.json"
    else
        echo "$SPECIES" >> "$LOG"
        rm -f "$TMP"
    fi

done < "$INPUT"