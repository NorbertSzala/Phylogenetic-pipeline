#!/usr/bin/env python3
#
# # Script to download proteome sequences from NCBI given on list with selected genomes from select_best_assemblies.py
# 3rd part of pipeline

# input: ranked selected_assemblies.tsv from select_best_assemblies.py
# Output: ../data/proteomes/sequences - .faa protein sequences


# ~~~~~ Imports ~~~~~
import subprocess
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
import argparse


# ~~~~~ Paths ~~~~~
parser = argparse.ArgumentParser(description="Download genomes basing on givel GCF id.")

parser.add_argument(
    "--input",
    type=Path,
    required=True,
    help="TSV file with ranked specie's RefSeqs assemblies",
)

parser.add_argument(
    "--output_zipped",
    required=True,
    help="Path where to save raw downloaded files from NCBI",
)

parser.add_argument(
    "--output_sequences",
    type=Path,
    required=True,
    help="Path where to save extracted .faa proteomes sequences",
)

args = parser.parse_args()
INPUT = Path(args.input)
OUTPUT_ZIPPED = Path(args.output_zipped)
OUTPUT_ZIPPED.mkdir(parents=True, exist_ok=True)

OUTPUT_SEQUENCES = Path(args.output_sequences)
OUTPUT_SEQUENCES.mkdir(parents=True, exist_ok=True)

LOG = OUTPUT_ZIPPED / "download_failed.log"
LOG.write_text("")  # clear log file


INCLUDE = "protein,gff3,genome,seq-report"


# ~~~~~ functions ~~~~~
def safe_name(name: str) -> str:
    """Change spaces and '/' to '__' - safe file names."""
    return name.replace(" ", "__").replace("/", "__")


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
    zip_path: Path, accession: str, species_safe: str, OUTPUT: Path
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
    out_faa = OUTPUT / f"{species_safe}.faa"
    out_faa.write_bytes(protein_path.read_bytes())

    return True


# ~~~~~ Read tsv input ~~~~~
print("#3. Downloading selected assemblies")

assemblies = defaultdict(list)


with INPUT.open() as fh:
    next(fh)  # leave header
    for line in fh:
        species, accession, rank = line.strip().split("\t")
        assemblies[species].append((int(rank), accession))

# ~~~~~ Downloading ~~~~~

for species, items in tqdm(
    assemblies.items(), total=len(assemblies), desc="Downloading assemblies"
):
    items.sort()  # rank 1 -> rank 2
    safe = safe_name(species)

    success = False

    for rank, accession in items:
        zip_path = OUTPUT_ZIPPED / f"{safe}_{accession}.zip"
        faa_path = OUTPUT_SEQUENCES / f"{safe}.faa"

        # if download successful, extract protein.faa
        if download(accession, zip_path):
            if extract_protein_faa(zip_path, accession, safe, OUTPUT_SEQUENCES):
                success = True
                break

            else:
                with LOG.open("a") as log:
                    log.write(f"NO_PROTEOME\t{species}\t{accession}\n")
        else:
            with LOG.open("a") as log:
                log.write(f"FAILED_RANK{rank}\t{species}\t{accession}\n")

    if not success:
        with LOG.open("a") as log:
            log.write(f"FAILED_ALL\t{species}\n")
