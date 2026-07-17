# PanTax-DBG: species- and strain-level profiling for paired-end short-read metagenomes

[![GitHub release](https://img.shields.io/github/v/release/xujialupaoli/PanTax-DBG)](https://github.com/xujialupaoli/PanTax-DBG/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pantax-dbg.svg)](https://anaconda.org/bioconda/pantax-dbg)
[![License](https://img.shields.io/github/license/xujialupaoli/PanTax-DBG)](LICENSE)

**PanTax-DBG** is a metagenomic profiler for estimating **species-level** and **strain-level** relative abundance from **paired-end short-read sequencing data**. It combines high-recall taxonomic candidate screening with colored compacted de Bruijn graph (ccDBG) refinement, enabling profiling across low- to high-abundance microbial signals.

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
  - [`my_genome_sizes.ncbi.tsv.gz`](#my_genome_sizesncbitsvgz)
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

### Install with Conda

PanTax-DBG 0.1.0 is available from Bioconda for Linux. Install it with `conda-forge` ahead of `bioconda` under strict channel priority:

```bash
conda create -n pantaxdbg_env pantax-dbg \
    --channel conda-forge \
    --channel bioconda \
    --strict-channel-priority

conda activate pantaxdbg_env

export PATH="$CONDA_PREFIX/libexec/pantax-dbg:$PATH"

pantax-dbg -h
```

The Conda package provides the PanTax-DBG command-line interface and all required runtime components. The commands above activate the environment and verify the installation.

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

Materialize the example genome paths for the current checkout:

```bash
awk -v root="$(pwd)" 'BEGIN {FS=OFS="\t"} NR==1 {print; next} {$5=root "/" $5; print}' \
    example/example_profile_genome_info.txt \
    > example/example_profile_genome_info.local.txt
```

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
    --ref-info example/example_profile_genome_info.local.txt \
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

The `--ref-info` table must correspond to the same genome collection and taxonomy used when building the database. Absolute genome paths are recommended for production analyses. The quick-start command above creates a local copy with absolute paths while leaving the distributed example table unchanged.

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
19703	1745244	1762272	1753758	2	custom_example
21691	1717340	2019581	1868460	2	custom_example
```

---

## Output files

The primary species and strain abundance tables report relative-abundance proportions. The hierarchical species tree retains the percentage scale of the taxonomic report, whereas the strain-refined tree reports proportions in its final column.

### `species_abundance.txt`

A tab-delimited species-level abundance table with a header.

| Column | Description |
| --- | --- |
| `speciesID` | Species-level taxonomic identifier. |
| `abundance` | Estimated relative abundance of the species. |

Example:

```text
speciesID	abundance
19703	0.5158330999999999
21691	0.4841669
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
9000023784	GENOME029482	0.2612231743
9000091645	GENOME127908	0.2546099257
9000024881	GENOME030794	0.2464849677
9000135004	GENOME204591	0.2376819323
```

This is the recommended file for downstream strain-level abundance analysis.

### `tax_profile.tre`

A headerless, tab-delimited hierarchical taxonomic abundance profile at non-strain ranks. It retains the nine-column cumulative report format.

| Column | Description |
| --- | --- |
| 1 | Taxonomic rank, such as `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, or `species`. |
| 2 | Taxonomic identifier. |
| 3 | Taxonomic lineage represented as taxonomic identifiers separated by pipe characters. |
| 4 | Taxon name. |
| 5 | Number of reads assigned uniquely to the taxon. |
| 6 | Number of reads assigned non-uniquely to the taxon after redistribution. |
| 7 | Number of assignments contributed by descendant taxa. |
| 8 | Cumulative assignments to the taxon and its descendants. |
| 9 | Cumulative relative abundance on a percentage scale. |

Example:

```text
root	1	1	root	0	0	240	240	100.00000
no rank	4	1|4	cellular organisms	0	0	240	240	100.00000
superkingdom	5	1|4|5	Archaea	0	0	120	120	51.58331
phylum	24	1|4|5|24	Thermoplasmatota	0	0	120	120	51.58331
```

### `tax_profile_strain.tre`

A headerless, tab-delimited hierarchical taxonomic profile after strain-level abundance refinement. This file includes strain-resolved branches and preserves a nine-column tree format.

| Column | Description |
| --- | --- |
| 1 | Taxonomic rank, including `strain` where applicable. |
| 2 | Taxonomic identifier. |
| 3 | Taxonomic lineage represented as taxonomic identifiers separated by pipe characters. |
| 4 | Taxon name. |
| 5 | Support/count field retained in the hierarchical tree representation. |
| 6 | Support/count field retained in the hierarchical tree representation. |
| 7 | Support/count field retained in the hierarchical tree representation. |
| 8 | Assigned or accumulated read-count field in the tree representation. |
| 9 | Final relative abundance after strain-level refinement. |

Example:

```text
root	1	1	root	0	0	240	240	1
no rank	4	1|4	cellular organisms	0	0	240	240	1
superkingdom	5	1|4|5	Archaea	0	0	120	120	0.5158331
phylum	24	1|4|5|24	Thermoplasmatota	0	0	120	120	0.5158331
class	217	1|4|5|24|217	Thermoplasmata	0	0	120	120	0.5158331
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

| Option | Default | Description |
| --- | --- | --- |
| `-r`, `--reads` | Required twice | Input paired-end FASTQ files, supplied in mate order. |
| `-d`, `--db-prefix` | Required | Prefix of a database generated by `pantax-dbg build-custom`. |
| `-i`, `--ref-info` | Required | Headered reference-genome information table for ccDBG refinement. |
| `-o`, `--output` | Required | Output directory for the final profiles. |
| `-t`, `--threads` | `8` | Number of worker threads. |
| `-k`, `--kmer-size` | `31` | k-mer size used during ccDBG refinement. |
| `--strain-min-reads` | `5.0` | Minimum assigned-read support before strain-candidate ranking. |
| `--strain-topk` | `5` | Maximum number of retained strain candidates per species. |
| `--strain-div` | `5.0` | Divisor used in secondary strain-support filtering. |
| `--r1-large-n` | `1000` | Candidate-set size above which large-set top-k filtering is applied. |
| `--r1-topk-large` | `3` | Number of candidates retained per species in large-set filtering mode. |
| `--r1-min-abundance` | `0.0` | Minimum first-round species abundance retained for ccDBG refinement. |
| `--species-min-abundance` | `0.0` | Minimum relative abundance retained in the final species profile. |

---

## Temporary files and local scratch space

During graph refinement, PanTax-DBG may create temporary intermediate files. Use `--ggcat-temp-dir` to direct temporary graph-building files to local scratch storage:

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

PanTax-DBG is currently distributed as versioned research software. If you use it in published work, please cite the software release:

```bibtex
@software{xu_pantax_dbg_2026,
  title   = {PanTax-DBG: species- and strain-level profiling for paired-end short-read metagenomes},
  author  = {Xu, Jialu},
  year    = {2026},
  month   = {6},
  version = {0.1.0},
  url     = {https://github.com/xujialupaoli/PanTax-DBG/releases/tag/v0.1.0},
  license = {GPL-3.0-or-later}
}
```

GitHub also provides APA and BibTeX formats through the **Cite this repository** link. A manuscript citation will be added after publication.

---

## License

PanTax-DBG is distributed under the GNU General Public License v3.0 or later. See [`LICENSE`](LICENSE) for details.
