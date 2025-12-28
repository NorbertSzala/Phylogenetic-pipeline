#!/usr/bin/env bash

# Script to select genome sequences from NCBI given a list of gene species.

# Returns the parameters for downloading genomes, like assembly, L50, N50 etc.
# 1st part of pipeline

#Input: list of species names in ../data/short_taxonomy.csv
# Output: ../data/proteomes/selection with selected assemblies for select the best assemblies in select_best_assemblies.py

# Next step: select_best_assemblies.py


# Safety settings - stop script on errors, undefined variables, or failed pipes
set -euo pipefail

source config.sh

# Set output and input directories/files
out_dir="../data/proteomes/selection"
mkdir -p "${out_dir}" ../logs
input_file="../data/short_taxonomy.csv"

# Set logging
LOG="../logs/selecting_assemblies.log"
: > "$LOG" # clear log file

# read species names line by line from the input file
while IFS= read -r SPECIES; do
    [[ -z "$SPECIES" ]] && continue

    SAFE=$(echo "$SPECIES" | tr ' /' '__')

    if ! datasets summary genome taxon "$SPECIES" \
        > "${out_dir}/${SAFE}_summary.json"; then
        echo "$SPECIES" | tee -a "$LOG"
        continue
    fi

done < "$input_file"