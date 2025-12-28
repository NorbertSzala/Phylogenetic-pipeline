#!/usr/bin/env python3
# 
# # Script to download proteome sequences from NCBI given on list with selected genomes from select_best_assemblies.py
# 3rd part of pipeline

# Input: ranked selected_assemblies.tsv from select_best_assemblies.py
# Output: ../data/proteomes/sequences - .faa protein sequences


import subprocess
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

# ~~~~~ Paths ~~~~~ 
INPUT_TSV = Path("../data/proteomes/selected_assemblies.tsv")
OUT_DIR = Path("../data/proteomes/metadata")
LOG_DIR = Path("../logs")
LOG_FILE = LOG_DIR / "download_failed.log"
SEQUENCE_DIR = Path("../data/proteomes/sequences")
SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)


INCLUDE = "protein,gff3,genome,seq-report"


OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.write_text("")  # clear log file

# ~~~~~ functions ~~~~~
def safe_name(name: str) -> str:
    """ Change spaces and '/' to '__' - safe file names. """
    return name.replace(" ", "__").replace("/", "__")



def download(accession: str, out_zip: Path) -> bool:
    """Download genome assembly from NCBI using datasets CLI.. Returns True/False."""
    # Create command
    cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--include", INCLUDE,
        "--filename", str(out_zip)
    ]
    result = subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL, #remove output, do not print to console
        stderr=subprocess.DEVNULL   #remove error output, do not print to console
    )
    return result.returncode == 0 # return true if download was successful


def extract_protein_faa(zip_path: Path, accession: str, species_safe: str, out_dir: Path) -> bool:
    """
    Unzip NCBI datasets archive and extract protein.faa.
    The file is renamed to <species_safe>.faa.
    """
    workdir = zip_path.parent

    # 1. Unzip archive
    result = subprocess.run(
        ["unzip", "-o", str(zip_path), "-d", str(workdir)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    if result.returncode != 0:
        return False

    protein_path = workdir / "ncbi_dataset" / "data" / accession / "protein.faa"

    if not protein_path.exists():
        return False

    # 2. Copy and rename protein.faa
    out_faa = out_dir / f"{species_safe}.faa"
    out_faa.write_bytes(protein_path.read_bytes())

    return True



# ~~~~~ Read tsv input ~~~~~
assemblies = defaultdict(list)


with INPUT_TSV.open() as fh:
    next(fh)  # leave header
    for line in fh:
        species, accession, rank = line.strip().split("\t")
        assemblies[species].append((int(rank), accession))

# ~~~~~ Downloading ~~~~~
for species, items in tqdm(assemblies.items(), total=len(assemblies), desc="Downloading assemblies"):
    items.sort()  # rank 1 -> rank 2
    safe = safe_name(species)

    success = False

    for rank, accession in items:
        zip_path = OUT_DIR / f"{safe}_{accession}.zip"
        faa_path = SEQUENCE_DIR / f'{safe}.faa'

        # if download successful, extract protein.faa
        if download(accession, zip_path):
            if extract_protein_faa(zip_path, accession, safe, SEQUENCE_DIR):
                success = True
                break
            
            else:
                with LOG_FILE.open("a") as log:
                    log.write(f"NO_PROTEOME\t{species}\t{accession}\n")
        else:
            with LOG_FILE.open("a") as log:
                log.write(f"FAILED_RANK{rank}\t{species}\t{accession}\n")

    if not success:
        with LOG_FILE.open("a") as log:
            log.write(f"FAILED_ALL\t{species}\n")