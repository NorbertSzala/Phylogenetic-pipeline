#!/usr/bin/env python3
"""
Compute Robinson–Foulds distance between two species trees.

- trees must have identical taxon labels
- branch lengths are ignored
"""

from ete3 import Tree
import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Compute RF distance between two trees")
    p.add_argument("--ref", type=Path, required=True, help="Reference tree")
    p.add_argument("--query", type=Path, required=True, help="Tree to compare")
    return p.parse_args()


def main():
    args = parse_args()

    ref = Tree(args.ref.read_text(), format=1)
    qry = Tree(args.query.read_text(), format=1)

    rf, max_rf, *_ = ref.robinson_foulds(qry, unrooted_trees=True)

    print("rf\tmax_rf\trf_norm")
    print(f"{rf}\t{max_rf}\t{rf / max_rf:.4f}")


if __name__ == "__main__":
    main()
