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
parser = argparse.ArgumentParser(
    description="Script creating mapping gene IDs to species name"
)

parser.add_argument(
    "--input",
    type=Path,
    required=True,
    help="Path to folder with downloaded in previous step .faa sequences",
)

parser.add_argument("--output", required=True, type=Path)

args = parser.parse_args()
INPUT = Path(args.input)
OUTPUT = Path(args.output)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
LOG = Path(OUTPUT.parent / "mapping_errors.log")
LOG.write_text("")


# ~~~~~ Functions ~~~~~


def create_mapping_gene_to_species(INPUT: Path, OUTPUT: Path):
    """Create mapping file with gene ID to species name."""
    if not INPUT.exists():
        raise SystemExit(f"Input directory does not exist: {INPUT}")

    fasta_files = list(INPUT.glob("*.faa"))
    if not fasta_files:
        raise SystemExit(f"No .faa files found in {INPUT}")

    with OUTPUT.open("w") as out:
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
                    with LOG.open("a") as log:
                        log.write(f"NO_HEADERS\t{fasta_file}\n")

            except Exception as e:
                with LOG.open("a") as log:
                    log.write(f"ERROR\t{fasta_file}\t{e}\n")


# ~~~~~ Main Logic ~~~~~
# Create mapping gene to species file
if __name__ == "__main__":
    create_mapping_gene_to_species(INPUT, OUTPUT)
