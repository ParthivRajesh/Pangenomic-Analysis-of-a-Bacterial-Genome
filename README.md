# Pangenomic Analysis of a Bacterial Genome (*Pantoea agglomerans*)

This project implements a complete pangenome analysis pipeline for Pantoea agglomerans, starting from raw NCBI genomes, annotating them with Prokka, and generating core-accessory gene clustering using Roary with downstream visualization. It provides an end-to-end, fully reproducible workflow.

## Tutorial Used

This project closely follows and expands on the pangenome tutorial:

- https://github.com/microgenomics/tutorials/blob/master/pangenome.md

---

## What is a Pangenome?

A **pangenome** refers to the full pool of genes found across all sequenced strains of a species. It includes:

- **Core genes** – present in (almost) all strains.
- **Soft core genes** – present in most strains, but not absolutely all.
- **Shell genes** – present in an intermediate fraction of strains.
- **Cloud genes** – rare genes that occur only in a few strains.

Together, these categories describe how conserved or variable the genome content is across a species.

---

## Organism: *Pantoea agglomerans*

- **Organism:** *Pantoea agglomerans*  
- **Type:** Gram-negative bacterium  
- **NCBI Taxon ID:** 549  
- **Reference article:** https://pmc.ncbi.nlm.nih.gov/articles/PMC1933083/ :contentReference[oaicite:5]{index=5}  

---

## Problem Statement

To annotate 52 genomes of *Pantoea agglomerans* and perform a pangenomic analysis using **ROARY** to identify:

- Core genes
- Soft core genes
- Shell genes
- Cloud genes

and interpret the pangenome structure of *Pantoea agglomerans* based on these categories. :contentReference[oaicite:6]{index=6}  

---

## Dataset

### NCBI Download

52 genomes of Pantoea agglomerans were downloaded manually from NCBI Datasets:

- NCBI Genome portal for P. agglomerans (Taxon ID 549). (https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=549)


#### Filters (before download)

- Annotated by NCBI RefSeq

- Assembly level: Chromosome and Complete

#### Download options (after clicking “Download”)

- Genome sequences (FASTA)

- Annotation features (GFF)

#### Metadata fields Extracted

- Organism scientific name
- Organism common name
- Organism qualifier
- Taxonomy ID
- Assembly name
- Assembly accession
- Source
- Annotation
- Level
- Contig N50
- Submission date
- Gene count
- BioProject
- BioSample 

---

## Observations

The pangenome structure of *Pantoea agglomerans* becomes clear when combining the numerical results from Roary with the visualizations generated using `roary_plots.py`.


### 1. Gene Frequency Distribution (Histogram)

The gene frequency histogram shows a strongly right-skewed pattern:

- Most genes occur in **only 1–5 genomes**, creating a dense cluster at the left side of the plot.  
- Only a very small number of genes are found in **all 52 genomes**, forming a tiny peak on the far right.

This pattern matches the summary statistics, where cloud genes dominate the pangenome (18,773 genes), while only 342 genes qualify as core. The histogram therefore highlights how large and variable the accessory genome is compared to the conserved component.


### 2. Accessory Gene Matrix + Phylogenetic Tree

The combined dendrogram and presence–absence matrix offer additional insight into strain relationships:

- The matrix is highly sparse, showing that most accessory genes are present in only a few strains.  
- Dark vertical bands represent groups of strains sharing certain gene clusters, but no widespread region of uniform presence exists outside the core genes.  
- The phylogenetic tree shows multiple clusters, indicating that strains with similar accessory gene profiles tend to group together, while the long branch lengths reflect high genomic diversity across the species.

Overall, the matrix and tree clearly demonstrate that *P. agglomerans* strains vary widely in their accessory gene content, consistent with the large cloud genome revealed by Roary.


### 3. Pangenome Composition (Pie Chart)

The pie chart provides an intuitive breakdown of gene categories:

- **Cloud genes (18,773)** occupy the overwhelming majority of the pangenome.  
- **Soft-core genes (2,773)** and **shell genes (1,813)** make up moderate proportions.  
- The **core genome (342 genes)** forms only a very small sliver of the total gene space.

This visual summary emphasizes the extremely large accessory genome and the strong genomic flexibility of the species.


### Observation Summary

Across all plots:

- The histogram shows most genes are rare, typical of an open pangenome.  
- The presence–absence matrix and tree highlight how much strain-to-strain variability is driven by accessory genes.  
- The pie chart visually confirms that cloud genes dominate the pangenome.  
- The Roary summary supports these findings with:
  - 342 core genes  
  - 2,773 soft-core genes  
  - 1,813 shell genes  
  - 18,773 cloud genes  
  - 23,701 total gene clusters  

Together, these observations show that *P. agglomerans* has **high genomic plasticity**, driven by a large collection of rare and niche-specific genes. This strongly supports the classification of its pangenome as **open**, with new genes likely to be discovered as more strains are sequenced.

---

## Results (ROARY) – Summary Statistics

After running Roary on the 52 *Pantoea agglomerans* genomes, the following gene category counts were obtained:

| Category         | Definition (Strain Presence)       | Gene Count |
|------------------|------------------------------------|-----------:|
| **Core genes**       | 99% ≤ strains ≤ 100%              | 342 |
| **Soft core genes**  | 95% ≤ strains < 99%               | 2,773 |
| **Shell genes**      | 15% ≤ strains < 95%               | 1,813 |
| **Cloud genes**      | 0% ≤ strains < 15%                | 18,773 |
| **Total genes**      | 0% ≤ strains ≤ 100%               | 23,701 |

These results indicate:

- A **small core genome** (342 genes) shared across almost all strains.  
- A **very large accessory (cloud) genome** (18,773 genes), indicating high genomic diversity and gene variability.  
- Significant numbers of **soft-core** and **shell** genes, representing moderately conserved gene families across multiple strains.

---

## Prerequisites

### 1. System & OS

- **OS:** Linux or macOS (tested on macOS).
- **CPU:** Multi-core (≥ 4 cores recommended).
- **RAM:** ≥ 8 GB (more is better for pangenomes).
- **Disk space:** A few GB for genomes, annotations, and Roary outputs.

On Apple Silicon, an x86-compatible conda environment (e.g. running under Rosetta) is recommended so that Prokka and Roary run smoothly. :contentReference[oaicite:7]{index=7}  

### 2. Software

- **conda** (or mamba) for environment management.
- **Prokka** for bacterial genome annotation.
- **Roary** for pangenome analysis.
- **Python 3** (compatible with `roary_plots.py`; Python 3 was used here).
- Python packages:
  - `matplotlib`
  - `seaborn`
  - `numpy`
  - `biopython`
- **wget** (or `curl`) for downloading `roary_plots.py`.
- Standard Unix tools: `bash`, `find`, `grep`, etc.

### 3. Suggested Conda Environment Setup

```bash
conda create -n bioinfo_x86 python=3.10 -y
conda activate bioinfo_x86

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install prokka roary -y
conda install matplotlib seaborn numpy biopython -y


