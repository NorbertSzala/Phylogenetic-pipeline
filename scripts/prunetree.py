#!/usr/bin/env python3
"""
Prune a phylogenetic tree to retain only taxa listed in a taxonomy file.

Input:
    - Tree file in Newick format
    - Text file with one species name per line

Output:
    - Pruned tree in Newick format

Matching rule:
    A leaf is kept if its name starts with the species name
    (spaces in species names are converted to underscores).

Design:
    - Pure tree manipulation (ETE)
    - No pipeline side effects
    - Deterministic and reproducible
"""

from ete4 import Tree
from pathlib import Path
import argparse
import sys


# Argument parsing


def parse_args():
    """
    Parse command-line arguments.

    Returns
    argparse.Namespace
        tree      : input tree file (Newick)
        taxonomy  : text file with species names
        output    : output pruned tree file
    """
    parser = argparse.ArgumentParser(
        description="Prune a phylogenetic tree using a list of species"
    )

    parser.add_argument(
        "--tree",
        type=Path,
        required=True,
        help="Input tree file in Newick format",
    )

    parser.add_argument(
        "--taxonomy",
        type=Path,
        required=True,
        help="Text file with species names (one per line)",
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output file for pruned tree (Newick)",
    )

    return parser.parse_args()


# Core functions


def load_tree(tree_file: Path) -> Tree:
    """
    Load a phylogenetic tree from file.

    Parameters
    tree_file : Path
        Path to Newick tree file

    Returns
    Tree
        ETE tree object
    """
    if not tree_file.exists():
        sys.exit(f"Tree file not found: {tree_file}")

    return Tree(str(tree_file))


def load_taxa(taxonomy_file: Path) -> list[str]:
    """
    Load species names from taxonomy file.

    Spaces are replaced with underscores to match tree labels.

    Parameters
    taxonomy_file : Path
        File with species names (one per line)

    Returns
    list[str]
        Normalized species names
    """
    if not taxonomy_file.exists():
        sys.exit(f"Taxonomy file not found: {taxonomy_file}")

    with taxonomy_file.open() as f:
        return [line.strip().replace(" ", "_") for line in f if line.strip()]


def find_matching_leaves(tree: Tree, taxa: list[str]) -> list[str]:
    """
    Identify leaf names that match the given taxa.

    A leaf is selected if its name starts with a taxon name.

    Parameters
    tree : Tree
        Input phylogenetic tree
    taxa : list[str]
        List of species identifiers

    Returns
    list[str]
        Names of leaves to keep
    """
    matched = []

    for node in tree.traverse():
        if node.is_leaf:
            for taxon in taxa:
                if node.name.startswith(taxon):
                    matched.append(node.name)
                    break

    return matched


def prune_tree(tree: Tree, leaves: list[str]) -> None:
    """
    Prune tree in place to retain selected leaves.

    Parameters
    tree : Tree
        Tree to prune
    leaves : list[str]
        Leaf names to keep
    """
    if not leaves:
        sys.exit("No matching leaves found — tree would be empty")

    tree.prune(leaves, preserve_branch_length=True)


def save_tree(tree: Tree, output_file: Path) -> None:
    """
    Save pruned tree to file.

    Parameters
    tree : Tree
        Pruned tree
    output_file : Path
        Output path
    """
    tree.write(outfile=str(output_file))


# Main


def main():
    args = parse_args()

    tree = load_tree(args.tree)
    taxa = load_taxa(args.taxonomy)

    matched_leaves = find_matching_leaves(tree, taxa)
    prune_tree(tree, matched_leaves)
    save_tree(tree, args.output)

    print(f"Pruned tree saved to: {args.output}")
    print(f"Leaves kept: {len(matched_leaves)}")


if __name__ == "__main__":
    main()
