# cbmf 🧬

**Combinatorial Bioinformatic Meta-Framework**:
automating everything from data acquisition to analysis.


## About

A collection of notes, scripts, and resources to automate bioinformatics workflows.

Code is designed to be modular and extensible.

[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)
Stop agonizing over stylistic details and just run prettier
over your markdown and html files.

```sh
prettier --write **/*.md **/*.html
```


### Structure of the Repository

Each folder contains a README.md with information about
the contents of the folder and how to use the scripts and resources.

The README.md files are written in GitHub Flavored Markdown (GFM), a superset of the standard Markdown syntax (note the `.md` extension), so they can be viewed either as plain text or as a rendered webpage on GitHub.

For more information, see the [GitHub Flavored Markdown Spec](https://github.github.com/gfm/). There is also a guide to [basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

which is a dialect of Markdown that is supported by GitHub. It extends the standard Markdown syntax with additional features that are useful for writing technical documentation. For more information, see the [GitHub Flavored Markdown Spec](https://github.github.com/gfm/). There is also a guide to [basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

This project uses GitHub Flavored Markdown (GFM) for documentation. GFM is a dialect of Markdown that is supported by GitHub that It extends the standard Markdown syntax with additional features that are useful for writing technical documentation. For more information, see the [GitHub Flavored Markdown Spec](https://github.github.com/gfm/). There is also a guide to [basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

> \[!TIP\]
> Use [markdownlint](https://github.com/DavidAnson/markdownlint) to lint markdown files.

## Workflow

### Data Acquisition

fastq files

#### Azenta

#### NextSeq550

#### MiSeq

### Data Processing

#### 0. Quality Control

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- another thing

#### 1. Alignment

- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

##### HISAT2

Copy, paste, and execute the following code to get started.

For more information, read the [manual](https://daehwankimlab.github.io/hisat2/manual/).

````sh

details, see the [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/).

```sh
#!/bin/bash

git clone "https://github.com/DaehwanKimLab/hisat2" && \
export PATH="$PATH:$(cd hisat2 && make -j \"$(nproc)\" && pwd)"
````

## Bookmarks

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

## Useful One-Liners

Update and upgrade everything on Ubuntu

```sh
sudo sh -c 'apt update && apt upgrade -y && apt dist-upgrade -y && apt autoremove -y && apt autoclean && apt clean'
```

Generate a txt file containing the md5sums of all files in a directory

```sh
find . -type f -exec md5sum {} \; > md5sums.txt
```

Check the files against the md5sums in the txt file

```sh
md5sum -c md5sums.txt
```

## `/genomes` Directory

The [`/genomes`](./genomes/README.md) directory contains the reference genomes and annotations for the organisms that are used in the pipelines. The files are downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) database.

Link to genomes \[link\]

## `samtools`

> Samtools is a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:
>
> - [Samtools](https:/github.com/samtools/samtools): Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
> - [BCFtools](https:/github.com/samtools/bcftools): Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants
> - [HTSlib](https:/github.com/samtools/): A C library for reading/writing high-throughput sequencing data

### key commands

- [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html) - sort alignments by leftmost coordinates

- [`samtools view`](https://www.htslib.org/doc/samtools-view.html) - converts between different formats

- [`samtools flagstat`](https://www.htslib.org/doc/samtools-flagstat.html) - quickly calculate simple statistics from a BAM file

- [`samtools index`](https://www.htslib.org/doc/samtools-index.html) - index a BAM file

- [`samtools merge`](https://www.htslib.org/doc/samtools-merge.html) - merge multiple sorted BAM files

- [`samtools mpileup`](https://www.htslib.org/doc/samtools-mpileup.html) - multi-way pileup

- [samtools](http://www.htslib.org/doc/samtools.html) is a suite of programs for interacting with high-throughput sequencing data.

- TODO: `REF_PATH` and `REF_CACHE`

### `.bam` to `.cram`

CRAM is a compressed version of the BAM format that is more efficient for long-term storage. It is a good idea to convert BAM files to CRAM files for long-term storage.

- [CRAM format specification](https://samtools.github.io/hts-specs/CRAMv3.pdf)
- [Using samtools to convert BAM to CRAM](https://www.htslib.org/workflow/cram.html)

**_IMPORTANT_**:

1. Alignments should be kept in chromosome/position sort order.
1. The reference must be available at all times. Losing it may be equivalent to losing all your read sequences.

when downloading the indexes for pipelines, see the compatibility issue below

```sh

$ samtools faidx download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz
[E::fai_build_core] File truncated at line 1
[E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
[faidx] Could not build fai index download/GCA_000001635.9_GRCm39_full_analysis_set.fna.gz.fai

```

decompress with gzip and recompress with bgzip

```sh
gunzip GCA_000001635.9_GRCm39_full_analysis_set.fna.gz && bgzip download/GCA_000001635.9_GRCm39_full_analysis_set.fna
```

## RNAseq

### Tuxedo Suite\[^1\]

1. HISAT2: A fast and sensitive alignment program for mapping next-generation sequencing reads (Kim et al., 2015)
1. StringTie: A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts (Pertea et al., 2015)
1. Ballgown: Flexible, isoform-level differential expression analysis (Frazee et al., 2015)

### Installation

Source code:

1. [HISAT2](htts://github.com/DaehwanKimLab/hisat2)
1. [StringTie](https://github.com/gpertea/stringtie)
1. [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html)

Function to install software and add to path:

```bash
install_and_add_to_path() {
  (
    git clone "$1" && cd "$(basename "$1" .git)" && make -j "$(nproc)"
    export PATH="$PATH:$(pwd)"
  )
}
```

> \[!CAUTION\]
> This function makes a lot of assumptions about the software being installed. It may not work for all software.

```bash
install_and_add_to_path https://github.com/DaehwanKimLab/hisat2.git
install_and_add_to_path https://github.com/gpertea/stringtie.git
```

### Ballgown

Start R and run:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ballgown")
```

### FastQC

Quality control of high throughput sequencing data.

#### Installation

```sh
sudo apt install openjdk-11-jdk &&
```

#### GTF (Gene Transfer Format)

- **Structure**: A specific version of GFF, often referred to as GFF2.
- **Attributes**: Contains nine fields per line, with the attributes field being a list of key-value pairs separated by semicolons.
- **Usage**: Primarily used in gene annotation pipelines and by databases like Ensembl.
- **Example**:

```gtf
chr1 HAVANA gene 11869 14409 . + . gene_id "ENSG00000223972"; gene_name "DDX11L1";
```

#### GFF (General Feature Format)

- **Versions**: Includes GFF2 and GFF3, with GFF3 being the more recent and widely used version.
- **Attributes**: Contains nine fields per line, with the attributes field formatted as a list of tag-value pairs separated by semicolons. GFF3 allows for hierarchical relationships between features.
- **Usage**: Commonly used for genome annotations and by databases like NCBI and UCSC Genome Browser.
- \*_Example_

```gff
##gff-version 3
chr1 . gene 11869 14409 . + . ID=gene00001;Name=DDX11L1;
```

- #### Key Differences
- **Format**: GTF is essentially a specific version of GFF (GFF2) with stricter formatting rules.
- **Versioning**: GFF3 supports more complex feature relationships and annotations compared to GTF/GFF2.
- **Usage**: Different databases and tools may prefer one format over the other based on their requirements and the complexity of the annotations.
  \*:
