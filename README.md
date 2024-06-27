# 🧬 Combinatorial Bioinformatic Meta-Framework

![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white)
![C++](https://img.shields.io/badge/c++-%2300599C.svg?style=for-the-badge&logo=c%2B%2B&logoColor=white)
![LaTeX](https://img.shields.io/badge/latex-%23008080.svg?style=for-the-badge&logo=latex&logoColor=white)
![Markdown](https://img.shields.io/badge/markdown-%23000000.svg?style=for-the-badge&logo=markdown&logoColor=white)
![Perl](https://img.shields.io/badge/perl-%2339457E.svg?style=for-the-badge&logo=perl&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

Notes, scripts, and resources to streamline and multiplex the analysis of
high-throughput sequencing data, from raw reads to biological insights.

## 📖 Table of Contents

1. 🔭 [Overview](#-overview)
2. 🚀 [Getting Started](#-getting-started)
3. 📚 [Documentation](#-documentation)
4. 📑 [Resources](#-resources)
5. [Workflows](#workflows)
   1. [Data Acquisition](#data-acquisition)
   2. [FastQC](#fastqc)
6. [One-liners](#useful-posix-compliant-one-liners)

## 🔭 Overview

This repository contains tools to automate key bioinformatic tasks:

- Data acquisition and storage
- Quality control
- Alignment to reference genomes
- Transcript assembly and quantification
- Differential expression analysis (work in progress)
- Visualization of results (work in progress)

## 🚀 Getting Started

If you are viewing these document on GitHub, you can copy the code snippets directly
by clicking the clipboard icon in the top right corner of code blocks like this one:

```sh
git clone https://github.com/rdnajac/cbmf.git
```

## 📚 Documentation

Check out the [documentation](./docs/README.md).

## 📑 Resources

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

## Workflows

This next section describes the steps involved in the analysis
of high-throughput sequencing data.

### Data Acquisition

#### Genomes

Stuff you should know:

- [FAQs - NCBI](https://ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/)
- GTF (a specific version of GFF2) and GFF (versions GFF2 and GFF3) are used in
  gene annotation, with GFF3 being more advanced.
- GTF uses key-value pairs for
  attributes, while GFF3 allows hierarchical feature relationships.

> [!INFO]
> The GenBank (GCA) assembly is an archival record that is owned by the submitter
> and may or may not include annotation.
> A RefSeq (GCF) genome assembly represents an NCBI-derived copy of a
> submitted GenBank (GCA) assembly.
> RefSeq (GCF) assembly records are maintained by NCBI.[^1]

##### NCBI Datasets and The Genome Reference Consortium

The most recent major releases from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets):

- [human](https://www.ncbi.nlm.nih.gov/grc/human)
- [mouse](https://www.ncbi.nlm.nih.gov/grc/mouse)

Link to the [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)

Download the data from the FTP server:

- [GRCm39 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/)
- [GRCh38 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)

##### `seqs_for_alignment_pipelines.ucsc_ids/`

> [!TIP]
> Skip building indexes from scratch by using the files
> provided in the `seqs_for_alignment_pipelines.ucsc_ids` folder.

| file                     | description            | action     |
| ------------------------ | ---------------------- | ---------- |
| fna.bowtie_index.tar.gz  | Bowtie2 index files    | `tar -xvf` |
| fna.hisat2_index.tar.gz  | HISAT2 index files     | `tar -xvf` |
| fna.fai                  | Samtools index file    | N/A        |
| fna.gz                   | FASTA format sequences | ?          |
| refseq_annotation.gff.gz | GFF3 format annotation | `gunzip`   |
| refseq_annotation.gtf.gz | GTF format annotation  | `gunzip`   |

#### Azenta

If you used Azenta for sequencing, they will send a link to directly download
the de-multiplexed FASTQ files (plus md5 checksums[^2]) from their sFTP server.
Click [here](https://3478602.fs1.hubspotusercontent-na1.net/hubfs/3478602/13012-M%26G%200222%20sFTP%20Guide-3.pdf)
for instructions on how to download data from Azenta's sFTP server.

Otherwise, consult the documentation for the appropriate Illumina sequencer:

- [MiSeq](https://support.illumina.com/sequencing/sequencing_instruments/miseq/documentation.html)
- [NextSeq500](https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550/documentation.html)

#### Align reads to a reference genome

##### Aligners

- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- STAR
- Subread

## Useful (POSIX-compliant) one-liners

Update and upgrade everything on Ubuntu

```sh
ys | sudo sh -c 'apt update && apt upgrade && apt dist-upgrade && apt autoremove && apt autoclean && apt clean'
```

Generate a txt file containing the md5sums of all files in a directory

```sh
md5sum ./*.fastq.gz > md5sums.txt
```

Generate a txt file containing the md5sums of all files in all subdirs

```sh
md5sum ./*/*.fastq.gz > checksums.md5
```

Check the files against the md5sums in the txt file

```sh
md5sum -c md5sums.txt
```

Export the current working directory to the PATH

```sh
export PATH=$PATH:$(pwd)
```

Check the total size of each folder in the current directory

```sh
du -sha --max-depth=1
```

The [`/genomes`](./genomes/README.md) directory contains information about
how to acquire and use reference genomes to align raw reads.

the reference genomes and annotations for the organisms that are used in the pipelines.
The files are downloaded from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly).

### [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[FastQC](https://github.com/s-andrews/FastQC) is a simple java application
that provides some quality control of high throughput sequencing data.

In order to run it, you need a suitable Java Runtime Environment (JRE) installed.

```sh
sudo apt install openjdk-11-jdk &&
```

Read the full installation instructions [here](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt).

[^1]: https://www.ncbi.nlm.nih.gov/genome/doc/assembly/
[^2]:
    The md5 checksum is a unique 32-character hexadecimal used to verify the
    integrity of the file during download or transfer. A checksum is computed for
    each file and changes if the file is modified. Read the original
    [RFC 1321](https://www.ietf.org/rfc/rfc1321.txt) for more information.
