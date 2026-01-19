#!/usr/bin/env python3
"""
Summarize filtering statistics and RF distance across pipeline runs.

python3 scripts/summarize_filtering_stats.py \
  --results results \
  --reference reference/reference_species_tree.tre \
  --output results/stats/summary

"""

from pathlib import Path
import pandas as pd
import argparse
from ete4 import Tree


# ---------------- CLI ----------------


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--results", type=Path, required=True, help="Root results directory")
    p.add_argument(
        "--reference", type=Path, required=True, help="Reference species tree"
    )
    p.add_argument("--output", type=Path, required=True)
    return p.parse_args()


# ---------------- RF ----------------


def rf_distance(ref_tree: Path, query_tree: Path) -> float:
    ref = Tree(ref_tree.read_text(), format=1)
    qry = Tree(query_tree.read_text(), format=1)
    rf, max_rf, *_ = ref.robinson_foulds(qry, unrooted_trees=True)
    return rf / max_rf if max_rf > 0 else None


# ---------------- Collect ----------------


def collect_runs(results_dir: Path, ref_tree: Path):
    records = []

    for run in results_dir.iterdir():
        if not run.is_dir():
            continue

        rec = {"run": run.name}

        # --- filtering stats ---
        stats = run / "stats/filtering"
        for mode in ["strict", "relaxed"]:
            f = stats / mode / "filter.stats.tsv"
            if not f.exists():
                continue

            df = pd.read_csv(f, sep="\t").set_index("metric")["value"]
            rec[f"{mode}_kept"] = int(df["kept_clusters"])
            rec[f"{mode}_input"] = int(df["input_clusters"])
            rec[f"{mode}_retention"] = rec[f"{mode}_kept"] / rec[f"{mode}_input"]

        # --- tree stats ---
        ts = run / "stats/tree_stats/tree_stats.tsv"
        if ts.exists():
            tdf = pd.read_csv(ts, sep="\t", names=["tree_id", "n_taxa"])
            rec["mean_n_taxa"] = tdf["n_taxa"].mean()
            rec["median_n_taxa"] = tdf["n_taxa"].median()
            rec["n_trees"] = len(tdf)

        # --- RF ---
        cons = run / "species_tree/consensus/consensus.treefile"
        if cons.exists():
            rec["rf_consensus"] = rf_distance(ref_tree, cons)

        astral = run / "species_tree/astral/astral.treefile"
        if astral.exists():
            rec["rf_astral"] = rf_distance(ref_tree, astral)

        records.append(rec)

    return pd.DataFrame(records)


def filter_bad_runs(df, n_genomes):
    """
    Remove biologically weak or unstable runs.
    """
    return df[
        (df["strict_kept"] >= 20)
        & (df["mean_n_taxa"] >= 0.8 * n_genomes)
        & (df["rf_consensus"] <= 0.3)
    ]


# ---------------- Main ----------------


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    df = collect_runs(args.results, args.reference)

    # remove weak runs
    df = filter_bad_runs(df, n_genomes=20)

    # rank remaining runs
    df = df.sort_values(
        by=["rf_consensus", "strict_kept", "mean_n_taxa"],
        ascending=[True, False, False],
    )

    df.to_csv(args.output / "ranked_runs.tsv", sep="\t", index=False)

    print("Summary written")
    print(df.head(10))


if __name__ == "__main__":
    main()
