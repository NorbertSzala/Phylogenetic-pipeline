#!/usr/bin/env python3
"""
Trim FASTA headers in alignment:
Species|GENE_ID  ->  Species

One alignment in, one alignment out.
"""

import argparse
from pathlib import Path
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    return p.parse_args()


def main():
    args = parse_args()

    if not args.input.exists():
        sys.exit(f"Input file not found: {args.input}")

    with args.input.open() as inp, args.output.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                header = line[1:].strip()
                species = header.split("|")[0]
                out.write(f">{species}\n")
            else:
                out.write(line)


if __name__ == "__main__":
    main()
