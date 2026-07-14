
# PanTax-DBG: species- and strain-level profiling for paired-end short-read metagenomes

[![GitHub release](https://img.shields.io/github/v/release/xujialupaoli/PanTax-DBG)](https://github.com/xujialupaoli/PanTax-DBG/releases)
[![License](https://img.shields.io/github/license/xujialupaoli/PanTax-DBG)](LICENSE)
<!-- Enable this badge after the Bioconda package is publicly available.
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pantax-dbg.svg)](https://anaconda.org/bioconda/pantax-dbg)
-->

**PanTax-DBG** is a metagenomic profiler for estimating **species-level** and **strain-level** relative abundance from **paired-end short-read sequencing data**. It combines high-recall taxonomic candidate screening with colored de Bruijn graph (DBG) refinement, enabling profiling across low- to high-abundance microbial signals.

> **Current public input type:** paired-end short reads (`-r read1.fq[.gz] -r read2.fq[.gz]`).

---

## Contents

- [Installation](#installation)
- [Commands](#commands)
- [Included example data](#included-example-data)
- [Quick start with the included example](#quick-start-with-the-included-example)
  - [1. Build a reference database](#1-build-a-reference-database)
  - [2. Profile paired-end reads](#2-profile-paired-end-reads)
- [Input file formats](#input-file-formats)
  - [`exampleDB_input_genomes.txt`](#exampledb_input_genomestxt)
  - [`example_profile_genome_info.txt`](#example_profile_genome_infotxt)
  - [`HRGM2_example.nodes.dmp` and `HRGM2_example.names.dmp`](#hrgm2_examplenodesdmp-and-hrgm2_examplenamesdmp)
  - [`my_genome_sizes.ncbi.tsv.gz`](#my_genome_sizesncbitasvgz)
- [Output files](#output-files)
  - [`species_abundance.txt`](#species_abundancetxt)
  - [`strain_abundance.txt`](#strain_abundancetxt)
  - [`tax_profile.tre`](#tax_profiletre)
  - [`tax_profile_strain.tre`](#tax_profile_straintre)
- [Common profiling options](#common-profiling-options)
- [Temporary files and local scratch space](#temporary-files-and-local-scratch-space)
- [Citation](#citation)
- [License](#license)

---

## Installation

### Install from Bioconda

Bioconda recommends using `conda-forge` with higher priority than `bioconda` and enabling strict channel priority. After the package is available in Bioconda:

```bash
conda create -n pantaxdbg_env pantax-dbg \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority

conda activate pantaxdbg_env

# Make the bundled PanTax-DBG backend tools visible.
export PATH="$CONDA_PREFIX/libexec/pantax-dbg:$PATH"

pantax-dbg -h
```

For general users, the packaged conda installation is recommended because it installs the required compiled components together with the command-line interface. A plain `pip install .` installation only installs the Python wrapper and is not sufficient unless the bundled native backends are also available through `PANTAX_DBG_LIBEXEC`.

If database construction stops because a backend executable cannot be found, the conda environment is active but its bundled backend directory is not on `PATH`. Run:

```bash
export PATH="$CONDA_PREFIX/libexec/pantax-dbg:$PATH"
```

and repeat the command. To make this permanent for the environment, the same line can be added to the user's shell startup file or to a conda activation script.

---

## Commands

PanTax-DBG provides two public commands:

```bash
pantax-dbg build-custom -h
pantax-dbg profile -h
```

| Command | Purpose |
| --- | --- |
| `pantax-dbg build-custom` | Build a custom taxonomic reference database from genome sequences and taxonomy files. |
| `pantax-dbg profile` | Estimate species- and strain-level relative abundance from paired-end short reads. |

---

## Included example data

The repository contains a compact runnable example under `example/`. The example is intentionally small (approximately 3 MB) and is meant for installation checks and command-line familiarization rather than biological benchmarking.

```text
example/
├── exampleDB_input_genomes.txt
├── example_profile_genome_info.txt
├── genome/
├── HRGM2_example.names.dmp
├── HRGM2_example.nodes.dmp
├── my_genome_sizes.ncbi.tsv.gz
├── README.md
├── read1.fq
└── read2.fq
```

| File or directory | Description |
| --- | --- |
| `genome/` | Four genome FASTA files used to build the toy example database. |
| `exampleDB_input_genomes.txt` | Headerless input manifest for `pantax-dbg build-custom`. |
| `HRGM2_example.nodes.dmp`, `HRGM2_example.names.dmp` | Taxonomy files used during database construction. |
| `my_genome_sizes.ncbi.tsv.gz` | Species genome-size information used for abundance calculation. |
| `read1.fq`, `read2.fq` | Small paired-end toy reads generated from the included genomes. |
| `example_profile_genome_info.txt` | Reference-genome information table used in profiling. |
| `README.md` | Step-by-step description of the included example. |
| `pantaxdbg_DB/` | Database directory created by the build command below. It is not required to be present before running the example. |
| `profile_res/` | Output directory created by the profiling command below. It is ignored by Git. |

---

## Quick start with the included example

Run the commands from the repository root:

```bash
cd PanTax-DBG

# Make the bundled PanTax-DBG backend tools visible.
export PATH="$CONDA_PREFIX/libexec/pantax-dbg:$PATH"
```

### 1. Build a reference database

```bash
mkdir -p example/pantaxdbg_DB

pantax-dbg build-custom \
    -k 31 \
    -w 51 \
    --verbose \
    --input-file example/exampleDB_input_genomes.txt \
    --taxonomy-files example/HRGM2_example.nodes.dmp example/HRGM2_example.names.dmp \
    --db-prefix example/pantaxdbg_DB/pantaxdbg_db \
    --level strain \
    --genome-size-files example/my_genome_sizes.ncbi.tsv.gz \
    -t 16 
```

The database prefix in this example is:

```text
example/pantaxdbg_DB/pantaxdbg_db
```

PanTax-DBG will generate database files associated with that prefix, including the files required for later profiling, such as:

```text
example/pantaxdbg_DB/pantaxdbg_db.hibf
example/pantaxdbg_DB/pantaxdbg_db.tax
```

For a large reference collection, database construction may require substantially more time and computational resources than this small tutorial.

### 2. Profile paired-end reads

```bash
mkdir -p example/profile_res

pantax-dbg profile \
    -r example/read1.fq \
    -r example/read2.fq \
    --db-prefix example/pantaxdbg_DB/pantaxdbg_db \
    --strain-min-reads 5.0 \
    --strain-topk 5 \
    --strain-div 5.0 \
    --r1-large-n 100 \
    --r1-topk-large 10 \
    --r1-min-abundance 1e-7 \
    --ref-info example/example_profile_genome_info.txt \
    --output example/profile_res \
    --threads 16 \
    -k 31 
```

The main outputs are:

```text
example/profile_res/species_abundance.txt
example/profile_res/strain_abundance.txt
example/profile_res/tax_profile.tre
example/profile_res/tax_profile_strain.tre
```

For a production run with more threads:

```bash
/usr/bin/time -v -o sample.time.log \
pantax-dbg profile \
    -r read1.fq.gz \
    -r read2.fq.gz \
    --db-prefix /path/to/pantaxdbg_db \
    --strain-min-reads 5.0 \
    --strain-topk 5 \
    --strain-div 5.0 \
    --r1-large-n 100 \
    --r1-topk-large 10 \
    --r1-min-abundance 1e-7 \
    --ref-info /path/to/example_profile_genome_info.txt \
    --output /path/to/sample_profile_res \
    --species-min-abundance 0.0 \
    --threads 64 \
    -k 31 \
    2>&1 | tee sample.profile.log
```

---

## Input file formats

### `exampleDB_input_genomes.txt`

This file is used by:

```bash
pantax-dbg build-custom --input-file exampleDB_input_genomes.txt
```

It is a **headerless, tab-delimited** file with three columns:

| Column | Name | Description |
| --- | --- | --- |
| 1 | `genome_path` | Path to a reference genome FASTA file. Files may be gzip compressed (`.fna.gz`). |
| 2 | `strain_name` | Unique name/identifier assigned to the genome or strain. |
| 3 | `strain_taxid` | Strain-level taxonomy identifier used in the database. |

Example:

```text
example/genome/GENOME000067.fna.gz	GENOME000067	9000000060
example/genome/GENOME000081.fna.gz	GENOME000081	9000000074
example/genome/GENOME000090.fna.gz	GENOME000090	9000000083
```

Requirements:

- One genome per line.
- `strain_name` should be unique.
- `strain_taxid` should be unique and consistent with the supplied taxonomy files.
- Genome paths must be accessible during database construction.

### `example_profile_genome_info.txt`

This file is used by:

```bash
pantax-dbg profile --ref-info example_profile_genome_info.txt
```

It is a **headered, tab-delimited** reference information table. The example file uses the following columns:

| Column | Required | Description |
| --- | --- | --- |
| `strain_name` | Yes | Name of the strain/genome. |
| `strain_taxid` | Yes | Strain-level taxonomic identifier; must match candidate strain identifiers in the database. |
| `species_taxid` | Yes | Species-level taxonomic identifier. |
| `species_name` | Recommended | Human-readable species name. |
| `genome_path` | Yes | Path to the FASTA sequence for the corresponding strain/genome. |

Example format:

```text
strain_name	strain_taxid	species_taxid	species_name	genome_path
GENOME000067	9000000060	42335	Species_A	example/genome/GENOME000067.fna.gz
GENOME000081	9000000074	43549	Species_B	example/genome/GENOME000081.fna.gz
```

The `--ref-info` table must correspond to the same genome collection and taxonomy used when building the database.

### `HRGM2_example.nodes.dmp` and `HRGM2_example.names.dmp`

These are taxonomy files supplied to:

```bash
--taxonomy-files HRGM2_example.nodes.dmp HRGM2_example.names.dmp
```

They provide the taxonomic structure and names for the database identifiers. The taxonomy files must contain entries compatible with the strain and species taxonomic identifiers used in the genome manifest and profiling reference table.

### `my_genome_sizes.ncbi.tsv.gz`

This gzip-compressed, tab-delimited file provides species-level genome-size statistics for abundance estimation and is supplied through:

```bash
--genome-size-files my_genome_sizes.ncbi.tsv.gz
```

Columns:

| Column | Description |
| --- | --- |
| `#species_taxid` | Species-level taxonomic identifier. |
| `min_ungapped_length` | Minimum ungapped genome length represented for that species. |
| `max_ungapped_length` | Maximum ungapped genome length represented for that species. |
| `expected_ungapped_length` | Expected ungapped genome length used for abundance normalization. |
| `number_of_genomes` | Number of genomes represented for that species. |
| `method_determined` | Method or source used to determine the genome-size record. |

Example:

```text
#species_taxid	min_ungapped_length	max_ungapped_length	expected_ungapped_length	number_of_genomes	method_determined
9000000530	2052496	2052496	2052496	1	custom
9000000531	2573894	2573894	2573894	1	custom
```

---

## Output files

All abundance values shown below are **relative-abundance proportions** rather than percentages. For example, `0.0601` represents approximately `6.01%` relative abundance.

### `species_abundance.txt`

A tab-delimited species-level abundance table with a header.

| Column | Description |
| --- | --- |
| `speciesID` | Species-level taxonomic identifier. |
| `abundance` | Estimated relative abundance of the species. |

Example:

```text
speciesID	abundance
42335	0.06012759102059925
43549	0.052470189404142255
```

This is the recommended file for downstream species-level abundance analysis.

### `strain_abundance.txt`

A tab-delimited strain-level abundance table with a header.

| Column | Description |
| --- | --- |
| `strain_taxid` | Strain-level taxonomic identifier. |
| `strain_name` | Genome/strain name in the reference information table. |
| `abundance` | Estimated relative abundance of the strain. |

Example:

```text
strain_taxid	strain_name	abundance
9000091645	GENOME127908	0.0001921999
9000023784	GENOME029482	0.0001633699
```

This is the recommended file for downstream strain-level abundance analysis.

### `tax_profile.tre`

A headerless, tab-delimited hierarchical taxonomic abundance profile at non-strain ranks.

| Column | Description |
| --- | --- |
| 1 | Taxonomic rank, such as `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, or `species`. |
| 2 | Taxonomic identifier. |
| 3 | Taxonomic lineage represented as taxonomic identifiers separated by `|`. |
| 4 | Taxon name. |
| 5 | Relative abundance of the taxon. |

Example:

```text
no rank	4	1|4	cellular organisms	0.9999999999999999
superkingdom	5	1|4|5	Archaea	0.0003555697272234092
phylum	24	1|4|5|24	Thermoplasmatota	0.0003555697272234092
class	217	1|4|5|24|217	Thermoplasmata	0.0003555697272234092
```

### `tax_profile_strain.tre`

A headerless, tab-delimited hierarchical taxonomic profile after strain-level abundance refinement. This file includes strain-resolved branches and preserves a nine-column tree format.

| Column | Description |
| --- | --- |
| 1 | Taxonomic rank, including `strain` where applicable. |
| 2 | Taxonomic identifier. |
| 3 | Taxonomic lineage represented as taxonomic identifiers separated by `|`. |
| 4 | Taxon name. |
| 5 | Support/count field retained in the hierarchical tree representation. |
| 6 | Support/count field retained in the hierarchical tree representation. |
| 7 | Support/count field retained in the hierarchical tree representation. |
| 8 | Assigned or accumulated read-count field in the tree representation. |
| 9 | Final relative abundance after strain-level refinement. |

Example:

```text
root	1	1	root	0	0	162083	162083	1.0000000011
no rank	4	1|4	cellular organisms	0	0	162083	162083	1.0000000011
superkingdom	5	1|4|5	Archaea	0	0	42	42	0.0003555698
phylum	24	1|4|5|24	Thermoplasmatota	0	0	42	42	0.0003555698
class	217	1|4|5|24|217	Thermoplasmata	0	0	42	42	0.0003555698
```

For most downstream statistical analyses, use `species_abundance.txt` or `strain_abundance.txt`; use the `.tre` files when a hierarchical taxonomic representation is required.

---

## Common profiling options

Run:

```bash
pantax-dbg profile -h
```

to view the complete command-line help available in your installed version.

Frequently used options:

| Option | Description |
| --- | --- |
| `-r` | Input read file. Supply exactly two paired-end read files by using `-r` twice. |
| `--db-prefix` | Prefix of a database previously generated by `pantax-dbg build-custom`. |
| `--ref-info` | Headered reference-genome information table for strain refinement. |
| `-o`, `--output`, `--out` | Output directory for profiling results. `--out` is kept as a compatibility alias for earlier scripts. |
| `--threads` | Number of threads used for computation. |
| `-k` | k-mer size used during graph-based refinement. |
| `--strain-min-reads` | Minimum read-support threshold applied during strain-level filtering. |
| `--strain-topk` | Maximum number of candidate strains retained per species in the corresponding refinement step. |
| `--strain-div` | Divisor used in secondary strain-support filtering. |
| `--r1-large-n` | Candidate-number threshold used to identify large candidate sets. |
| `--r1-topk-large` | Number of top candidates retained for large candidate sets. |
| `--r1-min-abundance` | Minimum abundance threshold used during first-round candidate selection. |
| `--species-min-abundance` | Minimum abundance threshold for final species-level reporting. |

---

## Temporary files and local scratch space

During graph refinement, PanTax-DBG may create temporary intermediate files. If your installed version exposes `--ggcat-temp-dir`, you can direct temporary graph-building files to local scratch storage:

```bash
pantax-dbg profile \
    -r read1.fq.gz \
    -r read2.fq.gz \
    --db-prefix /path/to/pantaxdbg_db \
    --ref-info /path/to/example_profile_genome_info.txt \
    --output /path/to/profile_res \
    --ggcat-temp-dir /tmp/pantax_dbg_${USER}_$$ \
    --threads 64 \
    -k 31
```

Using a local scratch filesystem such as `/tmp` can reduce I/O overhead on network-mounted storage. Make sure the selected location has sufficient available space for your sample and reference database.

---

## Citation

The citation for PanTax-DBG will be updated after publication of the associated manuscript. Until then, cite the software release and provide the repository URL:

```bibtex
@software{pantax_dbg_2026,
  title   = {PanTax-DBG: species- and strain-level profiling for paired-end short-read metagenomes},
  author  = {Xu, Jialu and Luo Group},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/xujialupaoli/PanTax-DBG}
}
```

---

## License

PanTax-DBG is distributed under the GNU General Public License v3.0 or later. See [`LICENSE`](LICENSE) for details.
