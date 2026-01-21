#!/usr/bin/env python3
# python strip_gene_ids_from_trees.py input.nwk output.species.nwk

from ete4 import Tree
import sys

inp = sys.argv[1]
out = sys.argv[2]

with open(inp) as fin, open(out, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue

        t = Tree(line)

        for leaf in t.leaves():
            leaf.name = leaf.name.split("|")[0]

        fout.write(t.write() + "\n")
