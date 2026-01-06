#!/usr/bin/env python3

"""
Filter weak gene trees based on bootstrap support.
One input tree -> keep or reject.
"""


# ~~~~~ Imports ~~~~~
from pathlib import Path
from collections import defaultdict
from statistics import mean, median
import argparse
import re
import sys


# ~~~~~ Arguments ~~~~~
parser = argparse.ArgumentParser(
    description="Filter a single gene tree based on bootstrap support"
)

parser.add_argument(
    "--input",
    type=Path,
    required=True,
    help="Input tree file (.treefile or .contree)",
)

parser.add_argument(
    "--output",
    type=Path,
    required=True,
    help="Output tree file (written only if tree is accepted)",
)

parser.add_argument(
    "--min_support",
    type=float,
    default=70,
    help="Minimal bootstrap support to consider node strong",
)

parser.add_argument(
    "--max_bad_frac",
    type=float,
    default=0.3,
    help="Max fraction of weak nodes allowed",
)


args = parser.parse_args()

INPUT = Path(args.input)
OUTPUT = Path(args.output)
MIN_SUPPORT = float(args.min_support)
MAX_BAD_FRAC = float(args.max_bad_frac)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
LOG = OUTPUT.parent / "rejecting_trees_failures.log"
LOG.write_text("")  # clear log file


# ~~~~~ Functions ~~~~~

# Bootstrap encoded like: )95:0.123 or )95.5:0.123
SUP_RE = re.compile(r"\)(\d+(?:\.\d+)?)(?=:)")


def supports_from_newick(s: str) -> list[float]:
    """Extract bootstrap supports from Newick string"""
    return [float(x) for x in SUP_RE.findall(s)]


def main():
    tree = INPUT.read_text().strip()
    supports = supports_from_newick(tree)

    # If no bootstrap values found -> reject (cannot evaluate)
    if not supports:
        sys.exit(1)

    bad = [x for x in supports if x < MIN_SUPPORT]  # list weak nodes in tree
    bad_frac = len(bad) / len(
        supports
    )  # evaluate fraction of weak nodes in trees. Avoid recejting tree with just one bad node

    # Reject weak tree
    if bad_frac > MAX_BAD_FRAC:
        sys.exit(1)

    # Tree accepted -> write output
    OUTPUT.write_text(tree + "\n")

    return 0


if __name__ == "__main__":
    main()
