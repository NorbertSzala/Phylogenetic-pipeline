#!/usr/bin/env python3

"""
Script to create a mapping file from gene IDs to species names
gene_001    Species_A

Usefull while analyzing clustering results from MMseqs2

#4th part of pipeline
"""

# ~~~~~ Imports ~~~~~
import argparse
from pathlib import Path
from tqdm import tqdm

# ~~~~~ Paths ~~~~~
INPUT_DIR = Path("../data/proteomes/sequences")
MAPPING_FILE = Path("../results/mapping/gene_to_species.tsv")
MAPPING_FILE.parent.mkdir(parents=True, exist_ok=True)

# ~~~~ Logging ~~~~~
LOG_FILE = Path("../logs/mapping_errors.log")
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
LOG_FILE.write_text("")


# ~~~~~ Functions ~~~~~


def create_mapping_gene_to_species(input_dir: Path, mapping_file: Path):
    """Create mapping file with gene ID to species name."""
    fasta_files = list(input_dir.glob("*.faa"))

    with mapping_file.open("w") as out:
        for fasta_file in tqdm(fasta_files, desc="Creating gene-to-species mapping"):
            species_name = fasta_file.stem
            found = False

            try:
                with fasta_file.open() as infile:
                    for line in infile:
                        if line.startswith(">"):
                            found = True
                            gene_id = line[1:].strip().split()[0]
                            out.write(f"{gene_id}\t{species_name}\n")

                if not found:
                    with LOG_FILE.open("a") as log:
                        log.write(f"NO_HEADERS\t{fasta_file}\n")

            except Exception as e:
                with LOG_FILE.open("a") as log:
                    log.write(f"ERROR\t{fasta_file}\t{e}\n")


# ~~~~~ Main Logic ~~~~~
# Create mapping gene to species file
if __name__ == "__main__":
    create_mapping_gene_to_species(INPUT_DIR, MAPPING_FILE)
