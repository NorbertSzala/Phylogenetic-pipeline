"""
Script to extract sequences for ONE orthologous gene cluster
and write them into a single FASTA file.

This script is designed to work with Nextflow:
- one cluster_id = one execution
- no alignments here (MAFFT is run by Nextflow)

Input:
- orthologs1to1.tsv
- all_proteomes.faa (indexed with faidx)

Output:
- one FASTA file with sequences for a single cluster
"""

# ~~~~~ Imports ~~~~~
import argparse
from pathlib import Path
from collections import defaultdict
from pyfaidx import Fasta
import sys


# ~~~~~ Argument parsing ~~~~~
parser = argparse.ArgumentParser(
    description="Extract sequences for a single orthologous cluster"
)

parser.add_argument(
    "--cluster_id",
    type=str,
    required=True,
    help="Cluster identifier (one cluster only)",
)

parser.add_argument(
    "--orthologs",
    type=Path,
    required=True,
    help="Path to orthologs1to1.tsv file",
)

parser.add_argument(
    "--proteomes",
    type=Path,
    required=True,
    help="Path to all_proteomes.faa (FASTA with all sequences)",
)

parser.add_argument(
    "--output",
    type=Path,
    required=True,
    help="Output FASTA file",
)

args = parser.parse_args()

CLUSTER_ID = args.cluster_id
ORTHOLOGS = args.orthologs
PROTEOMES = args.proteomes
OUTPUT = args.output

OUTPUT.parent.mkdir(parents=True, exist_ok=True)


# ~~~~~ Functions ~~~~~
def read_cluster_genes(tsv_file: Path, cluster_id: str) -> list[str]:
    """Eead orthologs1to1.tsv file and extract gene_ids
    belonging to a single cluster_id
    """
    genes = []

    with tsv_file.open() as f:
        next(f)  # skip header
        for line in f:
            cid, specie, gene_id = line.strip().split("\t")
            if cid == cluster_id:
                genes.append(gene_id)

    return genes


def write_cluster_fasta(
    cluster_id: str,
    genes: list[str],
    fasta_db: Fasta,
    output_fasta: Path,
):
    """Write sequences for a single cluster into FASTA file"""
    if len(genes) == 0:
        raise RuntimeError(f"No genes found for cluster {cluster_id}")

    if len(set(genes)) != len(genes):
        raise RuntimeError(f"Duplicate genes detected in cluster {cluster_id}")

    with output_fasta.open("w") as out:
        for gene_id in genes:
            try:
                seq = fasta_db[gene_id]
            except KeyError:
                raise RuntimeError(f"Sequence {gene_id} not found in FASTA file")

            out.write(f">{gene_id}\n")
            out.write(f"{str(seq)}\n")


# ~~~~~ Main logic ~~~~~

# open FASTA with its index
proteome_fasta = Fasta(
    PROTEOMES,
    as_raw=True,
    sequence_always_upper=True,
)

# extract gene list for this clustser
cluster_genes = read_cluster_genes(ORTHOLOGS, CLUSTER_ID)

# safety check (1-to-1 orthologs must have >1 sequence)
if len(cluster_genes) < 2:
    print(
        f"Cluster {CLUSTER_ID} has fewer than 2 sequences – skipping",
        file=sys.stderr,
    )
    sys.exit(1)

# write FASTA file
write_cluster_fasta(
    CLUSTER_ID,
    cluster_genes,
    proteome_fasta,
    OUTPUT,
)
