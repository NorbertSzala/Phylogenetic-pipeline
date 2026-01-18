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

| Tool       | Purpose                               | Tested version |
|-----------|----------------------------------------|----------------|
| Nextflow  | Workflow manager                       | ≥ 23.10       |
| Python    | Helper scripts                         | ≥ 3.9         |
| MMseqs2   | Sequence clustering                    | ≥ 14.7e284    |
| MAFFT     | Multiple sequence alignment            | ≥ 7.490       |
| IQ-TREE   | Gene trees & consensus species tree    | ≥ 2.2.0       |
| ASTRAL    | Species tree (summary method)          | 5.7.8         |
| Java      | Required for ASTRAL                    | ≥ 8           |

### Software availability

The pipeline assumes that all external tools are available in the user's `$PATH`,
except for ASTRAL, which is provided as a local JAR file.

Alternatively, absolute paths can be specified in `nextflow.config`.


## Installation (recommended: Conda)

Create a Conda environment with all dependencies:

```bash
conda create -n phylo \
  nextflow python=3.10 \
  mmseqs2 mafft iqtree openjdk \
  -c bioconda -c conda-forge
conda activate phylo
```

Install [Astral manually](https://github.com/smirarab/ASTRAL/)
```bash
wget https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
unzip Astral.5.7.8.zip
```

To test your installation, go to the place where you put the uncompressed ASTRAL, and run:
```bash
 java -jar astral.5.7.8.jar -i test_data/song_primates.424.gene.tre
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