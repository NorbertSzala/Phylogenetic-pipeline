#!/usr/bin/env python3

"""
Script to create a mapping file from gene IDs to species names
gene_001    Species_A

Usefull while analyzing clustering results from MMseqs2

parse one file per operation

#4th part of pipeline
"""

# ~~~~~ Imports ~~~~~
import argparse
from pathlib import Path
import sys


# ~~~~~ Paths ~~~~~
def parse_args():
    parser = argparse.ArgumentParser(
        description="Script creating mapping gene IDs to species name"
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to single .faa file",
    )

    parser.add_argument("--output", required=True, type=Path)

    return parser.parse_args()


# ~~~~~ Functions ~~~~~
def main():
    """Create mapping file with gene ID to species name."""

    args = parse_args()
    INPUT = args.input
    OUTPUT = args.output
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    if not INPUT.exists():
        raise SystemExit(f"Input file does not exist: {INPUT}")

    with OUTPUT.open("w") as out:
        species_name = INPUT.stem
        found = False

        try:
            with INPUT.open() as infile:
                for line in infile:
                    if line.startswith(">"):
                        found = True
                        gene_id = line[1:].strip().split()[0]
                        out.write(f"{gene_id}\t{species_name}\n")

            if not found:
                print(f'Current file {INPUT} does not contain ">".', file=sys.stderr)

        except Exception as e:
            sys.exit(f"Could not open {INPUT} file.")


# ~~~~~ Main Logic ~~~~~
# Create mapping gene to species file
if __name__ == "__main__":
    main()
