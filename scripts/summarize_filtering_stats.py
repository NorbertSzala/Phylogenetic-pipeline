#!/usr/bin/env python3

"""
Summarize ortholog filtering statistics across pipeline runs.

- Reads filter.stats.tsv files
- Aggregates strict / relaxed variants
- Produces:
  - summary table (CSV + TSV)
  - bar plots for report


Usage:
python3 scripts/summarize_filtering_stats.py \
  --input results/stats/filtering \
  --output results/stats/summary

"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import argparse


# --------------------------------------------------
# CLI
# --------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize filtering statistics from Nextflow pipeline"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/stats/filtering"),
        help="Root directory with filtering stats",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/stats/summary"),
        help="Output directory for summary tables and plots",
    )
    return parser.parse_args()


# --------------------------------------------------
# Main logic
# --------------------------------------------------


def collect_stats(root: Path) -> pd.DataFrame:
    records = []

    for stats_file in root.rglob("filter.stats.tsv"):
        df = pd.read_csv(stats_file, sep="\t")

        # infer metadata from path
        parts = stats_file.parts
        mode = "strict" if "strict" in parts else "relaxed"

        max_missing = None
        for p in parts:
            if p.startswith("max_missing_"):
                max_missing = int(p.replace("max_missing_", ""))

        df["mode"] = mode
        df["max_missing"] = max_missing
        df["path"] = str(stats_file)

        records.append(df)

    if not records:
        raise RuntimeError("No filter.stats.tsv files found")

    return pd.concat(records, ignore_index=True)


# --------------------------------------------------
# Plotting
# --------------------------------------------------


def plot_clusters(df: pd.DataFrame, outdir: Path):
    plt.figure(figsize=(7, 5))

    for mode in df["mode"].unique():
        sub = df[df["mode"] == mode]
        plt.bar(
            sub["max_missing"].fillna(0) + (0.2 if mode == "relaxed" else -0.2),
            sub["kept_clusters"],
            width=0.35,
            label=mode,
        )

    plt.xlabel("max_missing")
    plt.ylabel("Number of clusters")
    plt.title("Ortholog clusters retained after filtering")
    plt.legend()
    plt.tight_layout()

    plt.savefig(outdir / "clusters_retained.png", dpi=300)
    plt.close()


# --------------------------------------------------
# Entry point
# --------------------------------------------------


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    df = collect_stats(args.input)

    # Save tables
    df.to_csv(args.output / "filtering_summary.tsv", sep="\t", index=False)
    df.to_csv(args.output / "filtering_summary.csv", index=False)

    # Plots
    plot_clusters(df, args.output)

    print("✔ Filtering statistics summarized")
    print(f"  → {args.output / 'filtering_summary.tsv'}")
    print(f"  → {args.output / 'clusters_retained.png'}")


if __name__ == "__main__":
    main()
