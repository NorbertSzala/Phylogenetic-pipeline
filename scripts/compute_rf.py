#!/usr/bin/env python3
"""
Compute Robinson–Foulds distance between two species trees (ete4).

- trees must have identical taxon labels
- branch lengths are ignored
"""

from ete4 import Tree
import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Compute RF distance between two trees")
    p.add_argument("--ref", type=Path, required=True, help="Reference tree (Newick)")
    p.add_argument("--query", type=Path, required=True, help="Tree to compare (Newick)")
    return p.parse_args()


def main():
    args = parse_args()

    # read trees
    ref = Tree(args.ref.read_text())
    qry = Tree(args.query.read_text())

    # RF distance (unrooted)
    rf, max_rf, *_ = ref.robinson_foulds(qry, unrooted_trees=True)

    print("rf\tmax_rf\trf_norm")
    print(f"{rf}\t{max_rf}\t{rf / max_rf:.4f}")


if __name__ == "__main__":
    main()
