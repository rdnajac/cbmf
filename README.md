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

Installation instructions will generally be provided for Debian systems,
but the most scripts should work on any Unix-like system, e.g., macOS.

## Resources

## Structure of the Repository

Each folder contains a README.md with information about the folder's contents and
how to use the scripts and resources.
For more information about these documents, click [here](docs/README.md).

```plaintext
.
├── docs
│   ├── README.md

├── genomes
│   ├── README.md

├── scripts
│   ├── README.md

```

## Workflows

This next section describes the steps involved in the analysis
of high-throughput sequencing data.

### Data Acquisition

If you used Azenta for sequencing, they will send a link to directly download
the de-multiplexed FASTQ files (plus md5 checksums[^1]) from their sFTP server.
Click [here](https://3478602.fs1.hubspotusercontent-na1.net/hubfs/3478602/13012-M%26G%200222%20sFTP%20Guide-3.pdf)
for instructions on how to download data from Azenta's sFTP server.

Otherwise, consult the documentation for the appropriate Illumina sequencer:

- [MiSeq](https://support.illumina.com/sequencing/sequencing_instruments/miseq/documentation.html)
- [NextSeq500](https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550/documentation.html)

### bcl2fastq

> Read the [User Guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf).

If you have raw sequencing data in BCL format, you will need to convert it to
FASTQ format using the bcl2fastq2 Conversion Software.
This step can be skipped if you used Azenta for sequencing,
or if you correctly uploaded a valid sample sheet prior to sequencing.

#### Installation

Download bcl2fastq2 Conversion Software v2.20 Installer (Linux rpm) from
[Illumina](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html).

Check out this [post](https://www.biostars.org/p/266897/) for instructions
on how to convert this rpm (Red Hat Package Manager) file
into a deb (Debian Package Manager) file.

```sh
sudo alien -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
```

The `-i` flag installs the package after converting it to a temporary deb file.

### RNAseq

Tuxedo Suite\[^1\]

- [HISAT2](htts://github.com/DaehwanKimLab/hisat2):
  A fast and sensitive alignment program for mapping next-generation sequencing reads
- [StringTie](https://github.com/gpertea/stringtie):
  A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts
- [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html):
  Flexible, isoform-level differential expression analysis

\[^1\]: Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095

For more information, read the [manual](https://daehwankimlab.github.io/hisat2/manual/).

## Bookmarks

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

## Useful One-Liners

Update and upgrade everything on Ubuntu

```sh
ys | sudo sh -c 'apt update && apt upgrade && apt dist-upgrade && apt autoremove && apt autoclean && apt clean'
```

Generate a txt file containing the md5sums of all files in a directory

```sh
find . -type f -exec md5sum {} \; > md5sums.txt
```

Check the files against the md5sums in the txt file

```sh
md5sum -c md5sums.txt
```

Export the current working directory to the PATH

```sh
export PATH=$PATH:$(pwd)
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

[^1]:
    The md5 checksum is a unique 32-character hexadecimal used to verify the
    integrity of the file during download or transfer. A checksum is computed for
    each file and changes if the file is modified. Read the original
    [RFC 1321](https://www.ietf.org/rfc/rfc1321.txt) for more information.
