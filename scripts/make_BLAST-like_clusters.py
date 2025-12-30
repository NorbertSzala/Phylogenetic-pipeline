#!/usr/bin/env python3

# # Script to change names of genes in fasta headers and make BLAST, each sequence vs each sequence using MMseq2.
# Ive choosen MMseqs2 instead of Blast + MCL (markov cluster algorithm), becouse MMseqs does everything automatically and much faster
# 5th part of pipeline

# Input: protein sequences in .fasta format in data/proteomes/sequences
# Output:
# - merged and renamed fasta files (input for MMseqs2)
# - clusters_all__seqs.fasta unused
# - cluster_cluster.tsv - The most important file. Key to 3rd step - gene families. format: <member_sequence_id>    <representative_sequence_id>
# - clusters_rep_seqs.fasta - Representative sequences for clusters. one representative sequence for one cluster

# ~~~~~ Imports ~~~~~
import argparse
from pathlib import Path
from tqdm import tqdm
import subprocess
import shutil

# ~~~~~ Paths ~~~~~
INPUT_DIR = Path("../data/proteomes/sequences")
RESULTS_DIR = Path("../results/clusters/mmseqs2")
TMP_DIR = Path("../tmp/mmseqs")
LOG_FILE = Path("../logs/mmseqs2.log")

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
TMP_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
LOG_FILE.write_text("")


# ~~~~~ Functions ~~~~~


# At this moment it is unused, but can be useful in future
def rename_fasta(input_fasta: Path, output_fasta: Path):
    """Rename fasta headers to avoid conflicts.
    Add species name as prefix to each header."""
    species_name = input_fasta.stem

    with input_fasta.open() as infile, output_fasta.open("w") as outfile:
        for line in infile:
            if line.startswith(">"):
                header = line[1:].strip()
                new_header = f">{species_name}|{header}\n"
                outfile.write(new_header)

            else:
                outfile.write(line)


def merge_and_rename_fastas(input_dir: Path, merged_fasta: Path):
    """Merge and rename all fasta files in input_dir into a single fasta file."""
    fasta_files = list(input_dir.glob("*.faa"))
    if not fasta_files:
        raise RuntimeError("No FASTA files found in input directory")

    with merged_fasta.open("w") as outfile:
        for fasta in tqdm(fasta_files, desc="Merging fasta files", unit="files"):
            tmp_fasta = merged_fasta.parent / f"{fasta.stem}.tmp"
            rename_fasta(fasta, tmp_fasta)

            with tmp_fasta.open() as fh:
                outfile.write(fh.read())

            tmp_fasta.unlink()  # Remove temporary file


def run_mmseqs2(input_fasta: Path, tmp_dir: Path, output_dir: Path) -> bool:
    """RUn MMSeqs2 to create BLAST-like all-vs-all alignment database with clustering."""
    cmd = [
        "mmseqs",
        "easy-cluster",
        str(input_fasta),
        str(output_dir / "clusters"),
        str(tmp_dir),
        "--min-seq-id",
        "0.3",  # mimimum sequence identity 30%
        "-c",
        "0.8",  # coverage threshold 80% according to shorter sequence
        "--cov-mode",
        "1",
    ]

    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    with LOG_FILE.open("a") as log_file:
        log_file.write("STDOUT:\n")
        log_file.write(result.stdout)
        log_file.write("\nSTDERR:\n")
        log_file.write(result.stderr)

    return result.returncode == 0


# ~~~~~ Main logic ~~~~~
if __name__ == "__main__":

    merged_fasta = RESULTS_DIR / "all_proteomes.faa"

    merge_and_rename_fastas(INPUT_DIR, merged_fasta)

    print(f"Clustering sequences using MMseqs2")
    success = run_mmseqs2(merged_fasta, TMP_DIR, RESULTS_DIR)

    if not success:
        print("MMseqs2 failed. Check logs/mmseqs2.log")
        exit(1)

    # Cleanup only after successful run
    shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir()
