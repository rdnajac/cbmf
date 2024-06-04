# Combinatorial Bioinformatic Meta-Framework 🧬

Notes, scripts, and resources to streamline and multiplex the analysis of
high-throughput sequencing data, from raw reads to biological insights.

## Overview

This repository contains tools to automate key bioinformatic tasks:

- Data acquisition and storage
- Quality control
- Alignment to reference genomes
- Transcript assembly and quantification
- Differential expression analysis (work in progress)
- Visualization of results (work in progress)

These tools are in no way exhaustive,
but they provide a foundation for building more complex pipelines and workflows,
tailored to your experimental design and specific research questions.

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
the demultiplexed FASTQ files (plus md5 checksums[^1]) from their sFTP server.
Click [here](https://3478602.fs1.hubspotusercontent-na1.net/hubfs/3478602/13012-M%26G%200222%20sFTP%20Guide-3.pdf)
for instructions on how to download data from Azenta's sFTP server.

Othwerwise, consult the documentation for the appropriate Illumina sequencer:

- [MiSeq](https://support.illumina.com/sequencing/sequencing_instruments/miseq/documentation.html)
- [NextSeq500](https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550/documentation.html)

> [!NOTE]
> If you are using the BaseSpace command line interface (CLI),
> check out my quick reference guide [here](docs/basespace.md).

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

### Data Processing

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### HISAT2

Copy, paste, and execute the following code to get started.

For more information, read the [manual](https://daehwankimlab.github.io/hisat2/manual/).

```sh
git clone "https://github.com/DaehwanKimLab/hisat2" && \
export PATH="$PATH:$(cd hisat2 && make -j \"$(nproc)\" && pwd)"
```

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

The [`/genomes`](./genomes/README.md) directory contains information about
how to acquire and use reference genomes to align raw reads.

the reference genomes and annotations for the organisms that are used in the pipelines.
The files are downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) database.

## `samtools`

> Samtools is a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:
>
> - [Samtools](https:/github.com/samtools/samtools): Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
> - [BCFtools](https:/github.com/samtools/bcftools): Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants
> - [HTSlib](https:/github.com/samtools/): A C library for reading/writing high-throughput sequencing data

- [samtools](http://www.htslib.org/doc/samtools.html) is a suite of programs for interacting with high-throughput sequencing data.
- [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html) - sort alignments by leftmost coordinates
- [`samtools view`](https://www.htslib.org/doc/samtools-view.html) - converts between different formats
- [`samtools flagstat`](https://www.htslib.org/doc/samtools-flagstat.html) - quickly calculate simple statistics from a BAM file
- [`samtools index`](https://www.htslib.org/doc/samtools-index.html) - index a BAM file
- [`samtools merge`](https://www.htslib.org/doc/samtools-merge.html) - merge multiple sorted BAM files
- [`samtools mpileup`](https://www.htslib.org/doc/samtools-mpileup.html) - multi-way pileup

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

  <!-- Footnotes -->

[^1]:
    The md5 checksum is a unique 32-character hexadecimal used to verify the
    integrity of the file during download or transfer. A checksum is computed for
    each file and changes if the file is modified. Read the original
    [RFC 1321](https://www.ietf.org/rfc/rfc1321.txt) for more information.
