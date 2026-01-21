#!/usr/bin/env python3
# python select_core_taxon_set.py filtered25.nwk core25.nwk


from ete4 import Tree
from collections import Counter
import sys

inp = sys.argv[1]
out = sys.argv[2]

trees = []
sets = []

with open(inp) as f:
    for line in f:
        if not line.strip():
            continue
        t = Tree(line)
        s = frozenset(list(t.leaf_names()))
        trees.append(t)
        sets.append(s)

counter = Counter(sets)
core_set, core_count = counter.most_common(1)[0]

kept = 0
with open(out, "w") as fout:
    for t, s in zip(trees, sets):
        if s == core_set:
            fout.write(t.write() + "\n")
            kept += 1

print(f"Core taxa: {len(core_set)}", file=sys.stderr)
print(f"Trees in core: {core_count}", file=sys.stderr)
print(f"Written: {kept}", file=sys.stderr)
