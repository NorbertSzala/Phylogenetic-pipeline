#!/usr/bin/env bash
set -e

ASTRAL_VERSION=5.7.8
ASTRAL_DIR="tools/Astral"

mkdir -p ${ASTRAL_DIR}
cd ${ASTRAL_DIR}

if [ ! -f astral.${ASTRAL_VERSION}.jar ]; then
    echo "[INFO] Downloading ASTRAL ${ASTRAL_VERSION}"
    wget https://github.com/smirarab/ASTRAL/releases/download/v${ASTRAL_VERSION}/astral.${ASTRAL_VERSION}.jar
else
    echo "[INFO] ASTRAL already present"
fi
