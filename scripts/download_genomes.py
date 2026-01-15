#!/usr/bin/env python3
"""
Download proteome (.faa) from NCBI for one species based on ranked RefSeq assemblies.

Input:
    - One TSV file from select_best_assemblies.py
    (species, assembly_accession, rank)

Output:
    - One <species>.faa file written to output directory

Errors:
    - Reported via stderr
    - Non-zero exit code on failure (Nextflow-compatible)
"""

# ~~~~~ Imports ~~~~~
import subprocess
from pathlib import Path
from collections import defaultdict
import argparse
import sys
import shutil


# ~~~~~ Argument parsing ~~~~~
def parse_args():
    parser = argparse.ArgumentParser(
        description="Download proteome for one species from ranked RefSeq assemblies"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="TSV file with ranked RefSeq assemblies (one species)",
    )
    parser.add_argument(
        "--output_zipped",
        type=Path,
        required=True,
        help="Directory for downloaded NCBI zip files",
    )
    parser.add_argument(
        "--output_sequences",
        type=Path,
        required=True,
        help="Directory for output .faa proteome",
    )
    return parser.parse_args()


# ~~~~~ functions ~~~~~
INCLUDE = "protein,gff3,genome,seq-report"


def safe_name(name: str) -> str:
    """Change spaces and '/' to '__' - safe file names."""
    return name.replace(" ", "_").replace("/", "_")


def check_dependencies():
    if shutil.which("datasets") is None:
        sys.exit("Required program 'datasets' not found in PATH")
    if shutil.which("unzip") is None:
        sys.exit("Required program 'unzip' not found in PATH")


def download(accession: str, out_zip: Path) -> bool:
    """Download genome assembly from NCBI using datasets CLI.. Returns True/False."""
    # Create command
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        INCLUDE,
        "--filename",
        str(out_zip),
    ]
    result = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,  # remove output, do not print to console
        stderr=subprocess.DEVNULL,  # remove error output, do not print to console
    )
    return result.returncode == 0  # return true if download was successful


def extract_protein_faa(
    zip_path: Path, accession: str, species_safe: str, out_dir: Path
) -> bool:
    """
    Unzip NCBI datasets archive and extract protein.faa.
    The file is renamed to <species_safe>.faa.
    """
    workdir = zip_path.parent / accession
    workdir.mkdir(exist_ok=True)

    # 1. Unzip archive
    result = subprocess.run(
        ["unzip", "-o", str(zip_path), "-d", str(workdir)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    if result.returncode != 0:
        return False

    protein_path = workdir / "ncbi_dataset" / "data" / accession / "protein.faa"

    if not protein_path.exists():
        return False

    # 2. Copy and rename protein.faa
    out_faa = out_dir / f"{species_safe}.faa"
    with protein_path.open() as inp, out_faa.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                gene_id = line[1:].strip().split()[0]
                out.write(f">{species_safe}|{gene_id}\n")
            else:
                out.write(line)

    return True


# ~~~~~ Main logic ~~~~~


def main():
    args = parse_args()
    check_dependencies()  # check if all needed programmes are available

    INPUT = args.input
    OUT_ZIP = args.output_zipped
    OUT_SEQ = args.output_sequences

    OUT_ZIP.mkdir(parents=True, exist_ok=True)
    OUT_SEQ.mkdir(parents=True, exist_ok=True)

    assemblies = defaultdict(list)

    try:
        with INPUT.open() as fh:
            next(fh)  # leave header
            for line in fh:
                species, accession, rank = line.strip().split("\t")
                assemblies[species].append((int(rank), accession))

    except Exception as e:
        sys.exit(f"Failed to read TSV file {INPUT}: {e}")

    if not assemblies:
        sys.exit(f"No assemblies found in {INPUT}")

    # one species per TSV by design
    if len(assemblies) != 1:
        sys.exit(
            "Input TSV contains more than one species — "
            "this violates the pipeline contract"
        )

    species, items = next(iter(assemblies.items()))
    species_safe = safe_name(species)

    items.sort()  # rank 1 first

    for rank, accession in items:
        zip_path = OUT_ZIP / f"{species_safe}_{accession}.zip"

        if not download(accession, zip_path):
            print(
                f"Download failed (rank {rank}) for accession {accession}",
                file=sys.stderr,
            )
        if extract_protein_faa(zip_path, accession, species_safe, OUT_SEQ):
            # SUCCESS
            return

        print(
            f"No protein.faa found (rank {rank}) for accession {accession}",
            file=sys.stderr,
        )

    # if we reach h ere → all ranks failed
    sys.exit(f"All RefSeq assemblies failed for species: {species}")


if __name__ == "__main__":
    main()
