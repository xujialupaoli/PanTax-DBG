# PanTax-DBG

**PanTax-DBG** is a species- and strain-level metagenomic profiler that combines high-recall taxonomic pre-screening with colored de Bruijn graph based refinement.

PanTax-DBG bundles modified versions of ganon and ggcat, referred to here as **DBG-ganon** and **DBG-ggcat**. These modified backends are required by PanTax-DBG and are not identical to upstream ganon or ggcat. During installation, they are compiled from source and installed under `libexec/pantax-dbg/`. PanTax-DBG calls these bundled binaries directly rather than relying on system-wide ganon or ggcat installations.

<!-- [![BioConda Install](https://img.shields.io/conda/dn/bioconda/pantax_dbg.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/pantax_dbg)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pantax_dbg/badges/version.svg)](https://anaconda.org/bioconda/pantax_dbg)
[![License](https://img.shields.io/github/license/xujialupaoli/Themis)](https://www.gnu.org/licenses/gpl-3.0.en.html) -->


## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Features](#Features)
- [Quick start](#quick-start)
<!-- - [Input files](#input-files)
- [Output files](#output-files)
- [Command reference](#command-reference)
- [Tips and known issues](#tips-and-known-issues)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact) -->

---
## Overview
Themis is a fast and robust metagenomic profiler that achieves high accuracy across ultra-low to high sequencing depths. Themis combines a rapid, high-recall pre-screening step with graph-based refinement using colored de Bruijn graphs, reducing classification ambiguity and improving scalability to large reference databases.

## Installation 


```
conda create -n themis_env
conda activate themis_env
conda install -c bioconda -c conda-forge pantax_dbg
## Run pantax_dbg.
pantax_dbg -h
```
## Features

- **Commands**
  - `pantax-dbg build-custom` Build custom pantax_dbg databases.
  - `pantax-dbg profile` Profile reads against custom databases.


## Quick start
* **1-build-custom-reference-database** 
```
pantax-dbg build-custom  --input-file input_genomes.txt --taxonomy-files nodes.dmp names.dmp --db-prefix themisDB --level strain -t $threads -k 19 -w 51
```
input_genomes.txt is a headerless, tab-separated manifest where each line contains (1) the absolute path to a genome FASTA file, (2) its strain\_name, and (3) the corresponding strain-level NCBI taxid.

Due to the large size of the reference pangenome we used for testing, we provide the `genomes_info.txt` used here. You need to download these genomes from NCBI RefSeq and update the actual paths in `genomes_info.txt`. Please note that NCBI RefSeq periodically updates their database, so we cannot guarantee that all the listed genomes will be available. Building the reference pangenome takes approximately one week with this `genomes_info.txt`. 

* **2-profile**

```
# short read(pair-end)
pantax_dbg -r read1.fq -r $read2.fq --db-prefix themisDB --ref-info genomes_info.txt --out themis_query --threads 64 -k 31
# long read
pantax_dbg -r $reads.fq --single --db-prefix themisDB --ref-info genomes_info.txt --out themis_query --threads 64 -k 31
```
genomes_info.txt is a tab-separated metadata table with a header line. The columns are, in order: strain_name, strain_taxid, species_taxid, species_name, and genome_path, where strain_name and strain_taxid must be unique and genome_path gives the absolute path to the corresponding genome FASTA file.

Output file:species_abundance.txt
```
Taxonomic_ID    Relative_Abundance
12345           0.0012
...
```
Output file:tax_profile.tre
```
no rank 131567  1|131567        cellular organisms      1.0000000000000000
superkingdom    2       1|131567|2      Bacteria        1.0000000000000000
phylum  1224    1|131567|2|1224 Pseudomonadota  0.37078199931442135
class   1236    1|131567|2|1224|1236    Gammaproteobacteria     0.30406830971011906
order   135614  1|131567|2|1224|1236|135614     Xanthomonadales 0.006280138828970609
family  32033   1|131567|2|1224|1236|135614|32033       Xanthomonadaceae        0.006280138828970609
genus   68      1|131567|2|1224|1236|135614|32033|68    Lysobacter      0.006280138828970609
species 69      1|131567|2|1224|1236|135614|32033|68|69 Lysobacter enzymogenes  0.006280138828970609
...
```
