<!-- markdownlint-disable MD013 -->
# cbmf 🧬

**Combinatorial Bioinformatic Meta-Framework**: automating everything from data acquisition to analysis.

## About

This project is a collection of scripts and tools to automate bioinformatics workflows designed to be modular and extensible.

### Documentation

[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)

#### GitHub Flavored Markdown

This project uses GitHub Flavored Markdown (GFM) for documentation. GFM is a dialect of Markdown that is supported by GitHub that It extends the standard Markdown syntax with additional features that are useful for writing technical documentation. For more information, see the [GitHub Flavored Markdown Spec](https://github.github.com/gfm/). There is also a guide to [basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

> [!TIP]
> Use [markdownlint](https://github.com/DavidAnson/markdownlint) to lint markdown files.

## Bookmarks

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/) md
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

## Useful One-Liners

Update and upgrade everything on Ubuntu using `apt`:
```sh
sudo apt update && sudo apt upgrade -y && sudo apt dist-upgrade -y && sudo apt autoremove -y
```

## `/genomes` Directory

The [`/genomes`](./genomes/README.md) directory contains the reference genomes and annotations for the organisms that are used in the pipelines. The files are downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) database.

Link to genomes [link]

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
2. The reference must be available at all times. Losing it may be equivalent to losing all your read sequences.

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
