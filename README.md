# Phylogenomic pipeline (Nextflow DSL2)

Reproducible phylogenomic pipeline for inference of gene trees and species trees
using orthologous and pseudo-single-copy gene families.

## Overview
## Requirements
## Installation
## Running the pipeline
## Output structure
## Reproducibility notes
## Citation


``` bash
git clone https://github.com/NorbertSzala/Phylogenetic-pipeline
```


## Overview

This pipeline reconstructs gene trees and species trees starting from proteomes.

Main steps:
0. Create your species list.
1. Download proteomes from NCBI
2. All-vs-all sequence clustering (MMseqs2)
3. Ortholog filtering:
   - strict (1-to-1)
   - relaxed (N−k, pseudo-single-copy)
4. Multiple sequence alignment (MAFFT)
5. Gene tree inference (IQ-TREE)
6. Species tree inference:
   - consensus tree (IQ-TREE)
   - supertree / summary method (ASTRAL)


## Requirements

The pipeline requires the following software:

| Tool              | Purpose                             | Tested version       |
| ----------------- | ----------------------------------- | -------------------- |
| Nextflow          | Workflow manager                    | 25.10.2              |
| Python            | Runs scripts                        | 3.14.2               |
| ncbi-datasets-cli | Download data form NCBI             | 18.13.0              |
| MMseqs2           | Sequence clustering                 | 18-8cc5c             |
| MAFFT             | Multiple sequence alignment         | v7.525 (2024/Mar/13) |
| IQ-TREE           | Gene trees & consensus species tree | 2.0.7                |
| ASTRAL            | Species tree                        | 5.7.8                |
| Java              | Required for ASTRAL                 | 21.0.9               |


Software versions used in this study are recorded automatically in:

results/software/software_versions.txt


### Software availability

The pipeline assumes that all external tools are available in the user's `$PATH`,
except for ASTRAL, which is provided as a local JAR file.

Alternatively, absolute paths can be specified in `nextflow.config`.

## Scripts

The `scripts/` directory contains helper Python and Bash scripts used internally
by the pipeline (e.g. proteome download, ortholog filtering, FASTA extraction).
They are executed automatically by Nextflow and are not intended to be run manually.


## Installation (recommended: Conda)

Create a Conda environment with all dependencies and indall ASTRAL:

```bash
./scripts/install/install.sh

conda activate phylo_pipeline
```

## Input data

The pipeline expects a plain text file with one species name per line,
placed in:

data/short_taxonomy.csv

Example:

```bash
Saccharomyces cerevisiae
Malassezia globosa
Serpula lacrymans
Aspergillus niger
```

## Running the pipeline

Basic execution:

```bash
nextflow run main.nf
```

- Strict consensus species tree (1-to-1 orthologs):

``` bash
nextflow run main.nf \
  --consensus_mode strict
```

- Relaxed supertree using ASTRAL:

``` bash
nextflow run main.nf \
  --consensus_mode relaxed
```

- Enable bootstrap support for gene trees:

``` bash
nextflow run main.nf \
  --bootstrap true \
  --bootstrap_reps 1000
```


### Expected output structure:

```bash
tree ./
.
├── Astral                                  # Astral dependencies
├── data                                    # your file with species name should be placed there
├── results
│   ├── alignments                          # Alignments made by MMSeqs2
│   │   ├── relaxed
│   │   └── strict
│   ├── clusters                            # Clustering output. MMseqs2
│   │   ├── fasta
│   │   │   ├── relaxed
│   │   │   └── strict
│   │   └── mmseqs2
│   ├── gene_trees                          # Partial gene_trees - one tree per protein
│   │   ├── relaxed
│   │   └── strict
│   ├── mapping                             # Gene to species mapping
│   ├── proteomes                           # Proteome sequences downloaded form NCBI. Each organism have its own folder in ./sequences
│   │   └── sequences
│   └── species_tree                        # final ouptut
│       ├── astral
│       └── consensus
└── scripts                                 # Scripts to this pipeline

```

## Reproducibility notes

- The pipeline is deterministic except for ML tree inference.
- Random seeds are controlled internally by IQ-TREE when possible.
- All parameters used for a run are stored in Nextflow logs.
- The full pipeline can be resumed using `-resume`.


## Citation

If you use this pipeline, please cite:

- IQtree2:
    Bui Quang Minh, Heiko A Schmidt, Olga Chernomor, Dominik Schrempf, Michael D Woodhams, Arndt von Haeseler, Robert Lanfear, Corrigendum to: IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era, Molecular Biology and Evolution, Volume 37, Issue 8, August 2020, Page 2461, https://doi.org/10.1093/molbev/msaa131
- MMseqs2
    Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol. 2017 Nov;35(11):1026-1028. doi: 10.1038/nbt.3988. Epub 2017 Oct 16. PMID: 29035372.
- Astral
    Zhang, C., Rabiee, M., Sayyari, E. et al. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinformatics 19 (Suppl 6), 153 (2018). https://doi.org/10.1186/s12859-018-2129-y