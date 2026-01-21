#!/usr/bin/env python3
"""
Concatenate trimmed alignments (species-level) into a supermatrix.

All input FASTA files must contain sequences named by species.
Missing species are filled with gaps.
"""

import argparse
from pathlib import Path
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--species-order", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_fasta(path: Path) -> dict:
    seqs = {}
    current = None
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].strip()
                seqs[current] = ""
            else:
                seqs[current] += line.strip()
    return seqs


def load_species_order(path: Path) -> list[str]:
    with path.open() as f:
        return [line.strip() for line in f if line.strip()]


def main():
    args = parse_args()

    species_order = load_species_order(args.species_order)
    supermatrix = OrderedDict((sp, "") for sp in species_order)

    for aln_path in sorted(args.inputs):
        aln = read_fasta(aln_path)
        aln_len = len(next(iter(aln.values())))

        for sp in species_order:
            if sp in aln:
                supermatrix[sp] += aln[sp]
            else:
                supermatrix[sp] += "-" * aln_len

    with args.output.open("w") as out:
        for sp, seq in supermatrix.items():
            out.write(f">{sp}\n{seq}\n")


if __name__ == "__main__":
    main()
