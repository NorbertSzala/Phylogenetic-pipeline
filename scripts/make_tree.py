#!/usr/bin/env python3

"""
Script generating ML tree using IQtree.

input: single alignment of protein sequences

output: phylogenetic tree in .newick format

options: use bootstrap or not to compare results and reject the worst trees

"""


# ~~~~~ Imports ~~~~~
import argparse
from pathlib import Path
import subprocess

# ~~~~~ Paths ~~~~~
parser = argparse.ArgumentParser(
    description="Script counting trees by ML method using IQtree"
)

parser.add_argument(
    "--input",
    type=Path,
    required=True,
    help="Path to protein sequence alignment.",
)

parser.add_argument("--botstrap", required=True, type=int)
parser.add_argument("--output", required=True, type=Path)

args = parser.parse_args()
INPUT = Path(args.input)
BOOTSTRAP = int(args.bootstrap)
OUTPUT = Path(args.output)
OUTPUT.mkdir(parents=True, exist_ok=True)
LOG = Path(OUTPUT.parent / "mmseqs2_errors.log")
LOG.write_text("")
