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
from tqdm import tqdm

# ===== Paths =====
parser = argparse.ArgumentParser(
    description="Select best genome assemblies from NCBI RefSeq."
)
parser.add_argument(
    "--input",
    type=Path,
    required=True,
    help="Directory containing *_summary.json files",
)
parser.add_argument("--output", type=Path, required=True)


args = parser.parse_args()

INPUT = Path(args.input)
OUTPUT = Path(args.output)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
LOG = OUTPUT.parent / "NoRefSeq_assemblies.log"
LOG.write_text("")  # overwrite if exists


# ===== Helper functions =====


def assembly_priority(rep: dict) -> int:
    """
    Assign a simple priority based on assembly level.
    Higher value means better assembly.
    """
    level = rep.get("assembly_info", {}).get("assembly_level", "").lower()

    return {"chromosome": 3, "scaffold": 2, "contig": 1}.get(level, 0)


# ===== Main logic =====
print("#2. Selecting best assemblies")

with OUTPUT.open("w") as out:
    out.write("species\tassembly_accession\trank\n")

    for json_file in tqdm(
        INPUT.glob("*_summary.json"), desc="Selecting RefSeq assemblies"
    ):
        # Skip empty or invalid files
        if json_file.stat().st_size == 0:
            continue

        try:
            with json_file.open() as fh:
                data = json.load(fh)
        except json.JSONDecodeError:
            continue

        reports = data.get("reports", [])
        if not reports:
            continue

        # Species name is shared across reports in one file
        species = reports[0].get("organism", {}).get("organism_name", "UNKNOWN")

        # Keep only RefSeq assemblies (GCF)
        refseq = [
            rep for rep in reports if rep.get("current_accession", "").startswith("GCF")
        ]

        # If no RefSeq assembly exists, report species and skip it
        if not refseq:
            with LOG.open("a") as log:
                log.write(f"{species}\n")
            continue

        # Sort RefSeq assemblies by assembly level
        refseq.sort(key=assembly_priority, reverse=True)

        # Write top 2 RefSeq assemblies (main + backup)
        for rank, rep in enumerate(refseq[:2], start=1):
            acc = rep["current_accession"]
            out.write(f"{species}\t{acc}\t{rank}\n")
