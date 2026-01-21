#!/usr/bin/env python3
# python count_taxa_in_trees.py strict_gene_trees.species.nwk

from ete4 import Tree
from collections import Counter
import sys

inp = sys.argv[1]

counts = Counter()

with open(inp) as f:
    for line in f:
        if not line.strip():
            continue
        t = Tree(line)
        n = len(list(t.leaf_names()))
        counts[n] += 1

for k in sorted(counts):
    print(f"{k}\t{counts[k]}")
