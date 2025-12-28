#!/usr/bin/env python3

"""
Select RefSeq genome assemblies for downstream phylogenomic analysis.

Rationale:
Only assemblies with curated RefSeq annotations are used, as they
consistently provide annotated proteomes suitable for gene-based analyses.

Species without RefSeq assemblies are reported separately.
"""

import json
from pathlib import Path
from tqdm import tqdm

# ===== Paths =====

INDIR = Path("../data/proteomes/selection")
OUT_TSV = Path("../data/proteomes/selected_assemblies.tsv")
NO_REFSEQ = Path("../data/proteomes/no_refseq_species.txt")

NO_REFSEQ.write_text("")  # overwrite if exists


# ===== Helper functions =====

def assembly_priority(rep: dict) -> int:
    """
    Assign a simple priority based on assembly level.
    Higher value means better assembly.
    """
    level = rep.get("assembly_info", {}).get("assembly_level", "").lower()

    return {
        "chromosome": 3,
        "scaffold": 2,
        "contig": 1
    }.get(level, 0)


# ===== Main logic =====

with OUT_TSV.open("w") as out:
    out.write("species\tassembly_accession\trank\n")

    for json_file in tqdm(INDIR.glob("*_summary.json"), desc="Selecting RefSeq assemblies"):
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
        species = reports[0].get("organism", {}).get(
            "organism_name", "UNKNOWN"
        )

        # Keep only RefSeq assemblies (GCF)
        refseq = [
            rep for rep in reports
            if rep.get("current_accession", "").startswith("GCF")
        ]

        # If no RefSeq assembly exists, report species and skip it
        if not refseq:
            with NO_REFSEQ.open("a") as log:
                log.write(f"{species}\n")
            continue

        # Sort RefSeq assemblies by assembly level
        refseq.sort(key=assembly_priority, reverse=True)

        # Write top 2 RefSeq assemblies (main + backup)
        for rank, rep in enumerate(refseq[:2], start=1):
            acc = rep["current_accession"]
            out.write(f"{species}\t{acc}\t{rank}\n")