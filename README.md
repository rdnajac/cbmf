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

## About

[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)

This repository contains tools to automate key bioinformatic tasks:

- Data acquisition and storage
- Quality control
- Alignment to reference genomes
- Transcript assembly and quantification
- Differential expression analysis
- Visualization of results (work in progress)

There is also a wiki with additional resources and tutorials.
You can access the wiki by clicking on the tab at the top of the page
or by following [this link](https://github.com/rdnajac/cbmf/wiki).

## Usage

### Prerequisites

- Basic command line knowledge
- A POSIX-compliant shell (e.g. `bash`, `zsh`)
- A text editor (e.g. `vim`, `nano`, `emacs`)

```sh
git clone --recursive repo-url
```

For additional software, use micromamba for package management:

```sh
yes | "$SHELL" <(curl -L micro.mamba.pm/install.sh)
```

Piping `yes` into the script will accept the defaults:

```console
Micromamba binary folder? [~/.local/bin]
Init shell (bash)? [Y/n]
Configure conda-forge? [Y/n]
Prefix location? [~/micromamba]
```

Activate shell completion and restart the shell:

```sh
micromamba shell completion && exec "$SHELL"
```

### 💾 Data Acquisition

If you used a core facility, Azenta, or another commercial service for
sequencing, they will send link to directly download the de-multiplexed
FASTQ files, usually with corresponding md5 checksums.
For instructions on how to download data from Azenta's sFTP server,
click [here](https://3478602.fs1.hubspotusercontent-na1.net/hubfs/3478602/13012-M%26G%200222%20sFTP%20Guide-3.pdf).

> [!NOTE]
> sFTP (Secure File Transfer Protocol) provides an encrypted channel for data transfer.\
> The md5 checksum is a unique character sequence that is computed from the
> contents of a file and changes if the file is modified.
> Read the original [RFC 1321](https://www.ietf.org/rfc/rfc1321.txt).

Otherwise, consult the [documentation](https://developer.basespace.illumina.com/docs)
for the appropriate Illumina sequencer:

- [MiSeq](https://support.illumina.com/sequencing/sequencing_instruments/miseq/documentation.html)
- [NextSeq500](https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550/documentation.html)

#### BaseSpace and the `bs` command-line interface

The browser-based interface is useful for small-scale projects, but the
command-line interface is more efficient for large-scale projects.
Check out [examples](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples).

Install on macOS using Homebrew:

```sh
brew tap basespace/basespace && brew install bs-cli
```

After installation, authenticate with your BaseSpace credentials:

```sh
bs authenticate
```

Finally, run 'bs whoami' to verify that you are authenticated:

```text
+----------------+----------------------------------------------------+
| Name           | Ryan Najac                                         |
| Id             | ########                                           |
| Email          | rdn2108@cumc.columbia.edu                          |
| DateCreated    | 2021-07-13 15:29:51 +0000 UTC                      |
| DateLastActive | 2024-06-03 18:59:47 +0000 UTC                      |
| Host           | https://api.basespace.illumina.com                 |
| Scopes         | READ GLOBAL, CREATE GLOBAL, BROWSE GLOBAL,         |
|                | CREATE PROJECTS, CREATE RUNS, START APPLICATIONS,  |
|                | MOVETOTRASH GLOBAL, WRITE GLOBAL                   |
+----------------+----------------------------------------------------+
```

#### AWS S3: Simple Storage Service

s3 uris look like `s3://bucket-name/path/to/file.fastq.gz`

```sh
URI="s3://lab-aaf-ngs-data-archive/RNAseq/20240409_Tet2Rhoa-S1P1_RA/"
aws s3 sync "$URI" .
```

#### Demultiplexing Illumina sequencing data

Illumina instruments will demultiplex the data for you if you provide a valid
sample sheet prior to sequencing. For more information on how to create a sample
sheet, consult the [Illumina Support](https://support.illumina.com/) page.

Skip ahead to the [Quality Control](#quality-control) section if the
FASTQ files are already demultiplexed and ready for analysis, or keep reading
for instructions on how to demultiplex the data yourself.

Illumina hosts `.rpm` files for CentOS/RedHat Linux distros and the
source code (which must be compiled) for other distros.

Download bcl2fastq2 Conversion Software v2.20 Installer (Linux rpm) from
[Illumina](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html).

The AWS EC2 instance used for this project is based on Ubuntu, so we will
have to convert the `.rpm` file to a `.deb` file using the `alien` package,
as per this [post](https://www.biostars.org/p/266897/).

```sh
sudo alien -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
```

The following command is the default `bcl2fastq` command for demultiplexing
on the Nextseq, but with the `--no-lane-splitting` option added to combine
the reads from all four lanes into a single FASTQ file:

```sh
bcl2fastq --no-lane-splitting \
    --ignore-missing-bcls \
    --ignore-missing-filter \
    --ignore-missing-positions \
    --ignore-missing-controls \
    --auto-set-to-zero-barcode-mismatches \
    --find-adapters-with-sliding-window \
    --adapter-stringency 0.9 \
    --mask-short-adapter-reads 35 \
    --minimum-trimmed-read-length 35 \
    -R "$run_folder" \
    -o "$output_folder" \
    --sample-sheet "$sample_sheet" \
```

> [!WARNING]
> As of (DATE?), `bcl2fastq` is no longer supported; use `bclconvert` instead.\
> You can install `bclconvert` using the same methods as described above.

Read the docs:

- [bcl2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf)
- [bclconvert](https://support-docs.illumina.com/SW/BCL_Convert_v4.0/Content/SW/BCLConvert/BCLConvert.htm)

#### Quality Control

| Tool                                                                 | Description                                               | Installation                                                                                                           |
| -------------------------------------------------------------------- | --------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Generates html reports containing straightforward metrics | [source](https://github.com/s-andrews/FastQC)                                                                          |
| [GATK](https://gatk.broadinstitute.org/hc/en-us)                     | Analyzes high-throughput sequencing data                  | [instructions](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt)                                  |
| [Picard Tools](https://broadinstitute.github.io/picard/)             | Manipulates high-throughput sequencing data               | comes packaged with [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) |

To run these QC applications, you need a suitable Java Runtime Environment (JRE):

```sh
sudo apt-get install -y openjdk-11-jdk
```

> [!TIP]
> After aligning the reads to the reference genome, these tools can be re-ran on the
> resulting BAM files to ensure that the alignment was successful or to consolidate
> the results from paired-end sequencing.

#### Reference Genomes

The reference genome is a digital nucleotide sequence that represents the
genome of an organism. It is used as a reference for mapping reads from
high-throughput sequencing experiments.

The most recent major releases from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets):

- [human](https://www.ncbi.nlm.nih.gov/grc/human)
- [mouse](https://www.ncbi.nlm.nih.gov/grc/mouse)

Link to the [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)

Download the data from the FTP server:

- [GRCm39 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/)
- [GRCh38 (latest major release) FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/)

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

- [FAQs - NCBI](https://ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/)
- GTF (a specific version of GFF2) and GFF (versions GFF2 and GFF3) are used in
  gene annotation, with GFF3 being more advanced.
- GTF uses key-value pairs for
  attributes, while GFF3 allows hierarchical feature relationships.
- The GenBank (GCA) assembly is an archival record that is owned by the submitter
  - may or may not include annotation.
- A RefSeq (GCF) genome assembly represents an NCBI-derived copy of a submitted GenBank (GCA) assembly.
  - RefSeq (GCF) assembly records are maintained by [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/assembly/)

### 📑 Additional Resources

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

## 🧬 Align reads to a reference genome

> In this study, the algorithmically different mappers...
> were used to map experimentally generated RNA-Seq data from the...
> higher plant _Arabidopsis thaliana_ and to quantify the transcripts.[^1]

[^1]: Schaarschmidt et al., 2020

| Tool     | Description                                                                              |
| -------- | ---------------------------------------------------------------------------------------- |
| Bwa      | Maps short DNA sequences to reference genome; Uses BWT for indexing                      |
| STAR     | Specialized for RNA-Seq; uses seed-extension search; detects splice-junctions.           |
| HISAT2   | Splice-aware; uses graph-based alignment for DNA and RNA sequences.                      |
| RSEM     | Quantifies transcript abundances; uses expectation-maximization and pre-defined mappers. |
| salmon   | Uses quasi-mapping; BWT, suffix array and FMD algorithm for shared substring discovery.  |
| kallisto | Uses pseudo-alignments; maps k-mers to De Bruijn graph for isoform quantification.       |
| CLC      | Read mapping approach by Mortazavi et al.; only commercial tool with a GUI.              |

Other Aligners:

- [HISAT2](https://daehwankimlab.github.io/hisat2/)
  - splice-aware and suitable for transcriptome based or RNA-seq alignment
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - suitable for genome based alignment (ChIP-seq/WGS)
- [Subread](https://subread.sourceforge.net/)[^4], also comes packaged with:
  - [Subjunc](https://subread.sourceforge.net/subjunc.html)[^5], an exon-exon junction detector
  - [featureCounts](https://subread.sourceforge.net/featureCounts.html), a read summarization program

[^4]:
    The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote,
    **_Y Liao, GK Smyth, W Shi_**, Nucleic acids research, 2013 [PMID:23558742](https://pubmed.ncbi.nlm.nih.gov/23558742/)

[^5]:
    featureCounts: an efficient general purpose program for assigning sequence reads to genomic features,
    **_Y Liao, GK Smyth, W Shi_**, Bioinformatics, 2014 [PMID:24227677](https://pubmed.ncbi.nlm.nih.gov/24227677/)

## RNAseq

The transcriptome is the complete set of RNA transcripts in a cell,
and RNA-seq is used to analyze the transcriptome by sequencing the RNA molecules
after reverse transcription to cDNA.

### Tuxedo Suite

Transcript-level expression analysis of RNA-seq experiments
with HISAT, StringTie and Ballgown.[^2]

[^2]: Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. <https://doi.org/10.1038/nprot.2016.095>

- [Protocol](https://www.nature.com/articles/nprot.2016.095)
- [Software Availability](https://ccb.jhu.edu/software.shtml)

| Tool      | Description                                                | Manual                                                                      | Source                                             |
| --------- | ---------------------------------------------------------- | --------------------------------------------------------------------------- | -------------------------------------------------- |
| HISAT2    | Align reads to reference genome                            | [manual](https://daehwankimlab.github.io/hisat2/manual/)                    | [source](htts://github.com/DaehwanKimLab/hisat2)   |
| StringTie | Assembler of RNA-Seq alignments into potential transcripts | [manual](https://ccb.jhu.edu/software/stringtie/index.shtml)                | [source](https://github.com/gpertea/stringtie)     |
| Ballgown  | Flexible, isoform-level differential expression analysis   | [manual](https://bioconductor.org/packages/release/bioc/html/ballgown.html) | [source](https://github.com/alyssafrazee/ballgown) |

> [!TIP]
> StringTie come packaged with `gffcompare` for comparing and evaluating the
> accuracy of RNA-seq transcript assemblers. Read the
> [manual](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) for details.

#### Workflow

The following steps correspond to the workflow described in the protocol:

- Steps 1-2: Align reads to the reference genome and sort to BAM format
- Steps 3-6: `stringtie`
- Steps 7-9: `ballgown` setup
- Steps 10-15: Differential expression analysis
  - filter to remove low-abundance genes
  - identify transcripts
  - identify genes
  - add gene names
  - sort by p-value
  - write to file
- Steps 16-18: Visualization

Tips:

- Save the bg object to a file for later use: `save(bg, file="bg.rda")`
- Add log2 fold change and p-value to the gene-level data frame:
  `log2fc <- log2(gene$mean2/gene$mean1)`

## ChIPseq

Chromatic immunoprecipitation sequencing (ChIP-seq) is a method used to analyze
protein interactions with DNA. It combines chromatin immunoprecipitation (ChIP)
with massively parallel DNA sequencing to identify the binding sites of
DNA-associated proteins.

### DROMPA

## Acknowledgements

Shout out to these awesome docs:

- [Learn Vimscript the Hard Way](https://learnvimscriptthehardway.stevelosh.com/)
- [tao-of-tmux](https://tao-of-tmux.readthedocs.io/)
- [mamba](https://mamba.readthedocs.io/)

<!-- References -->

<!-- Bioconda -->

[^1]: Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences. Nature Methods, 2018 doi:10.1038/s41592-018-0046-7.

<!-- conda-env-mod -->

[^2]: A. K. Maji, L. Gorenstein and G. Lentner, "Demystifying Python Package Installation with conda-env-mod," 2020 IEEE/ACM International Workshop on HPC User Support Tools (HUST) and Workshop on Programming and Performance Visualization Tools (ProTools), GA, USA, 2020, pp. 27-37, doi: 10.1109/HUSTProtools51951.2020.00011.
