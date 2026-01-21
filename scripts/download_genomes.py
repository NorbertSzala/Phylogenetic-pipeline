#!/usr/bin/env python3
"""
Download and validate proteomes (.faa) from NCBI for one species based on ranked RefSeq assemblies.

Input:
    - One TSV file from select_best_assemblies.py
    (species, assembly_accession, rank)

Output:
    - One <species>.faa file written to output directory

Errors:
    - Reported via stderr
    - Non-zero exit code on failure (Nextflow-compatible)

Validation steps:
1. Robust downloading with fallback over ranked assemblies
2. Extraction of protein FASTA
3. FAST proteome-level quality control (QC)
4. Transparent logging of all decisions


Design principles:
- stdout  : structured, human-readable progress
- stderr  : reasons for rejection / errors
- outputs : ONLY data products (FASTA + TSV), never logs
- exit 0  : proteome accepted
- exit 1  : all assemblies rejected

"""

# ~~~~~ Imports ~~~~~
import subprocess
from pathlib import Path
from collections import defaultdict
import argparse
import sys
import shutil


# ~~~~~ Paremeters ~~~~~
MIN_PROTEINS = 4000
MAX_PROTEINS = 30000
MIN_MEAN_LEN = 250
MIN_PROTEIN_LEN = 30
MAX_SHORT_FRAC = 0.25

NCBI_INCLUDE = "protein"


# ~~~~~ Argument parsing ~~~~~
def parse_args():
    """
    Parse command-line arguments.

    Expected input:
    - TSV with one species and ranked RefSeq assemblies
    - output directories for archives and FASTA files
    - path to log file (important for reproducibility)
    """
    parser = argparse.ArgumentParser(
        description="Download proteome and make quality control for one species from ranked RefSeq assemblies"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="TSV file with ranked RefSeq assemblies (one species)",
    )
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument("--min_proteins", type=int, default=4000)
    parser.add_argument("--max_proteins", type=int, default=30000)
    parser.add_argument("--min_mean_len", type=int, default=250)
    parser.add_argument("--min_protein_len", type=int, default=30)
    parser.add_argument("--max_short_frac", type=float, default=0.25)

    return parser.parse_args()


# ~~~~~ functions ~~~~~
INCLUDE = "protein,gff3,genome,seq-report"


def safe_species_name(name: str) -> str:
    """Change spaces and '/' to '__' - safe file names."""
    return name.replace(" ", "_").replace("/", "_")


def check_dependencies():
    """Check required programmes"""
    for prog in ["datasets", "unzip"]:
        if shutil.which(prog) is None:
            sys.exit(f"Required program not found: {prog}")


def run(cmd: list[str]) -> bool:
    return (
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,  # remove output, do not print to console
            stderr=subprocess.DEVNULL,  # remove error output, do not print to console
        ).returncode
        == 0
    )  # return true if download was successful


def download_assembly(accession: str, out_zip: Path) -> bool:
    """Download genome assembly from NCBI using datasets CLI.. Returns True/False."""
    # Create command
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        NCBI_INCLUDE,
        "--filename",
        str(out_zip),
    ]
    return run(cmd)


def extract_protein_faa(zip_path: Path, accession: str) -> Path | None:
    """
    Unzip NCBI datasets archive and extract protein.faa.
    The file is renamed to <species_safe>.faa.
    """
    workdir = zip_path.parent / accession
    workdir.mkdir(exist_ok=True)

    # 1. Unzip archive
    if not run(["unzip", "-o", str(zip_path), "-d", str(workdir)]):
        return None

    faa = workdir / "ncbi_dataset" / "data" / accession / "protein.faa"
    return faa if faa.exists() else None


def protein_lengths(faa: Path) -> list[int]:
    """Read FASTA and return protein lengths only."""
    lengths = []
    seq = ""

    with faa.open() as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    lengths.append(len(seq))
                seq = ""
            else:
                seq += line.strip()
        if seq:
            lengths.append(len(seq))

    return lengths


def validate_fungal_proteome(faa: Path) -> tuple[bool, str]:
    """
    Fast biological QC for fungal proteomes.

    Returns
    -------
    (accepted, reason)
    """
    lengths = protein_lengths(faa)
    n = len(lengths)

    if n < MIN_PROTEINS:
        return False, f"too few proteins ({n})"
    if n > MAX_PROTEINS:
        return False, f"too many proteins ({n})"

    mean_len = sum(lengths) / n
    if mean_len < MIN_MEAN_LEN:
        return False, f"mean length too small ({mean_len:.1f})"

    short_frac = sum(l < MIN_PROTEIN_LEN for l in lengths) / n
    if short_frac > MAX_SHORT_FRAC:
        return False, f"too many short proteins ({short_frac:.1%})"

    return True, "OK"


# ~~~~~ Main logic ~~~~~


def main():
    args = parse_args()
    check_dependencies()

    args.outdir.mkdir(parents=True, exist_ok=True)

    qc_tsv = args.outdir / "qc_summary.tsv"

    # ------------------------------------------------------------------
    # Load assemblies (one species per TSV by pipeline contract)
    # ------------------------------------------------------------------
    assemblies = defaultdict(list)
    with args.input.open() as f:
        next(f)
        for line in f:
            species, acc, rank = line.strip().split("\t")
            assemblies[species].append((int(rank), acc))

    if len(assemblies) != 1:
        print("No valid assemblies for this species — skipping", file=sys.stderr)
        return 0

    species, ranked = next(iter(assemblies.items()))
    species_safe = safe_species_name(species)
    ranked.sort()

    # Try assemblies in rank order
    with qc_tsv.open("w") as qc:
        qc.write("species\taccession\trank\tstatus\treason\n")

        for rank, acc in ranked:
            print(f"[{species}] trying {acc} (rank {rank})")

            zip_path = args.outdir / f"{species_safe}_{acc}.zip"

            if not download_assembly(acc, zip_path):
                print(f"[{species}] download failed: {acc}", file=sys.stderr)
                qc.write(f"{species}\t{acc}\t{rank}\trejected\tdownload_failed\n")
                continue

            faa = extract_protein_faa(zip_path, acc)
            if faa is None:
                print(f"[{species}] protein.faa missing: {acc}", file=sys.stderr)
                qc.write(f"{species}\t{acc}\t{rank}\trejected\tno_protein_faa\n")
                continue

            ok, reason = validate_fungal_proteome(faa)
            if not ok:
                print(f"[{species}] rejected {acc}: {reason}", file=sys.stderr)
                qc.write(f"{species}\t{acc}\t{rank}\trejected\t{reason}\n")
                continue

            # ACCEPTED
            out_faa = args.outdir / f"{species_safe}.faa"
            with faa.open() as inp, out_faa.open("w") as out:
                for line in inp:
                    if line.startswith(">"):
                        gid = line[1:].split()[0]
                        out.write(f">{species_safe}|{gid}\n")
                    else:
                        out.write(line)

            qc.write(f"{species}\t{acc}\t{rank}\taccepted\tOK\n")
            print(f"[{species}] accepted {acc}")
            return 0

    # All assemblies failed
    print(f"[{species}] all assemblies rejected", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
