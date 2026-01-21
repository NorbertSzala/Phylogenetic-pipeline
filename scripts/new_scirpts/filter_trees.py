#!/usr/bin/env python3
# python filter_trees_by_taxa_count.py input.nwk filtered25.nwk 25

from ete4 import Tree
import sys

inp = sys.argv[1]
out = sys.argv[2]
TARGET = int(sys.argv[3])

kept = 0

with open(inp) as fin, open(out, "w") as fout:
    for line in fin:
        if not line.strip():
            continue
        t = Tree(line)
        if len(list(t.leaf_names())) == TARGET:
            fout.write(t.write() + "\n")
            kept += 1

print(f"Kept trees: {kept}", file=sys.stderr)
