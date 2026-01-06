#!/usr/bin/env python3

"""
Select RefSeq genome assemblies for downstream phylogenomic analysis.

Only assemblies with RefSeq annotations are used (accession starting with GCF).

Species without RefSeq assemblies are reported separately.

2nd part of pipeline
"""

# ~~~~~ Imports ~~~~~
import argparse
import json
from pathlib import Path
import sys

# ===== Paths =====
parser = argparse.ArgumentParser(
    description="Select best RefSeq assemblies from one NCBI summary JSON."
)

parser.add_argument(
    "--input",
    type=Path,
    required=True,
    help="Single *_summary.json file",
)

parser.add_argument("--output", type=Path, required=True, help="Output TSV file")


args = parser.parse_args()

INPUT = Path(args.input)
OUTPUT = Path(args.output)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)


# ===== Helper functions =====
def assembly_priority(rep: dict) -> int:
    """
    Assign a simple priority based on assembly level.
    Higher value means better assembly.
    """
    level = rep.get("assembly_info", {}).get("assembly_level", "").lower()

    return {"chromosome": 3, "scaffold": 2, "contig": 1}.get(level, 0)


def load_summary(path: Path) -> dict:
    """Function returning loaded JSON fil"""
    if not path.exists():
        sys.exit(f"Input file does not exists: {path}")
    if path.stat().st_size == 0:
        sys.exit(f"Empty input file: {path}")

    try:
        with path.open() as fh:
            return json.load(fh)
    except json.JSONDecodeError:
        sys.exit(f"Invalid JSON: {path}")


def select_refseq_assemblies(data: dict) -> tuple[str, list[dict]]:
    """Funciton selecting the best refseq assemblies to download later"""
    reports = data.get("reports", [])
    if not reports:
        return "UNKNOWN", []

    species = reports[0].get("organism", {}).get("organism_name", "UNKNOWN")

    refseq = [
        rep for rep in reports if rep.get("current_accession", "").startswith("GCF")
    ]

    refseq.sort(key=assembly_priority, reverse=True)
    return species, refseq


def write_output(
    output: Path,
    species: str,
    refseq: list[dict],
    max_hits: int = 2,
) -> None:
    with output.open("w") as out:
        out.write("species\tassembly_accession\trank\n")
        for rank, rep in enumerate(refseq[:max_hits], start=1):
            acc = rep["current_accession"]
            out.write(f"{species}\t{acc}\t{rank}\n")


def main() -> None:
    data = load_summary(INPUT)
    species, refseq = select_refseq_assemblies(data)
    write_output(OUTPUT, species, refseq)


if __name__ == "__main__":
    main()
