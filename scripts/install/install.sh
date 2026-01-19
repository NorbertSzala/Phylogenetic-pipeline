#!/usr/bin/env bash
set -e

ENV_NAME=phylo_pipeline

echo "[INFO] Creating conda environment: ${ENV_NAME}"

if command -v mamba &> /dev/null; then
    mamba env create -f install/environment.yml || mamba env update -f install/environment.yml
else
    conda env create -f install/environment.yml || conda env update -f install/environment.yml
fi

echo "[INFO] Activating environment"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${ENV_NAME}

echo "[INFO] Installing ASTRAL"
bash install//install_astral.sh

echo "[INFO] Installation finished"
echo "Run:"
echo "  conda activate ${ENV_NAME}"
# echo "  nextflow run main.nf"
