python3 - <<'PY'
from ete3 import Tree

inp = "manual_treebuild/strict_gene_trees.nwk"
out = "manual_treebuild/strict_gene_trees.species.nwk"

with open(inp) as fin, open(out, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue

        t = Tree(line, format=1)

        for leaf in t:
            # usuń wszystko po |
            leaf.name = leaf.name.split("|")[0]

        fout.write(t.write(format=1) + "\n")

print("Zapisano drzewa z nazwami gatunków:")
print(out)
PY
