# 🧬 Combinatorial Bioinformatic Meta-Framework

Notes, scripts, and resources to streamline and multiplex the analysis of
high-throughput sequencing data, from raw reads to biological insights.

## 📚 Table of Contents

- [🔭 Overview](#-overview)
- [🚀 Getting Started](#-getting-started)
  - [📦 Package Management](#-package-management)
- [🧬 Sequencing Workflows](#-sequencing-workflows)
  - [💾 Data Acquisition](#-data-acquisition)
  - [🔀 Demultiplexing](#-demultiplexing)
  - [🔍 Quality Control](#-quality-control)
  - [🏗️ Alignment](#-alignment)
  - [🔬 Assembly and Quantification](#-assembly-and-quantification)
  - [📊 Visualization](#-visualization)
  - [📈 Differential Expression](#-differential-expression)
- [📑 Additional Resources](#-additional-resources)

## 🔭 Overview

This repository contains tools to automate key bioinformatic tasks:

- Data acquisition and storage
- Quality control
- Alignment to reference genomes
- Transcript assembly and quantification
- Differential expression analysis
- Visualization of results
- Package management and dependency resolution

There is also a wiki with additional resources and tutorials.
You can access the wiki by clicking on the tab at the top of the page
or by following [this link](https://github.com/rdnajac/cbmf/wiki).

## 🚀 Getting Started

Clone the repository:

```sh
git clone https://github.com/rdnajac/cbmf.git
```

### 📦 Package Management

The only other thing we need to install is `micromamba`, a lightweight package manager
that can be used to install `conda` environments, like bioconda[^1].

[^1]: Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences. Nature Methods, 2018 doi:10.1038/s41592-018-0046-7.

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

## 🧬 Sequencing Workflows

The following sections outline the steps involved in analyzing high-throughput sequencing data.

A quick note on file formats:

| File Type              | Description                                 | Typical Extension                    | Details                                                        |
| ---------------------- | ------------------------------------------- | ------------------------------------ | -------------------------------------------------------------- |
| FASTQ                  | Raw reads from sequencer                    | .fastq, .fastq.gz                    | Human readable; often compressed (.gz)                         |
| FASTA                  | Sequence data                               | .fasta, .fa, .fasta.gz, .fa.gz       | Human readable; often compressed (.gz)                         |
| GTF                    | Gene Transfer Format (GFF2 variant)         | .gtf, .gtf.gz                        | Human readable; can be compressed                              |
| GFF2                   | General Feature Format v2                   | .gff, .gff2, .gff.gz, .gff2.gz       | Human readable; can be compressed                              |
| GFF3                   | General Feature Format v3                   | .gff3, .gff3.gz                      | Human readable; can be compressed                              |
| SAM                    | Sequence Alignment/Map                      | .sam, .sam.gz                        | Human readable; can be compressed                              |
| BAM                    | Binary Alignment/Map                        | .bam                                 | Binary (compressed); not human readable                        |
| CRAM                   | Compressed Reference-oriented Alignment Map | .cram                                | Binary (compressed); not human readable                        |
| VCF                    | Variant Call Format                         | .vcf, .vcf.gz                        | Human readable; often compressed                               |
| BED                    | Browser Extensible Data                     | .bed, .bed.gz                        | Human readable; can be compressed                              |
| TSV                    | Tab-Separated Values                        | .tsv, .tsv.gz                        | Human readable; can be compressed                              |
| CSV                    | Comma-Separated Values                      | .csv, .csv.gz                        | Human readable; can be compressed                              |
| TXT                    | Plain Text                                  | .txt, .txt.gz                        | Human readable; can be compressed                              |
| Bowtie2 Index          | Bowtie2 aligner index                       | .bt2, .bt2.tar.gz                    | Binary; not human readable; often compressed                   |
| BWA Index              | BWA aligner index                           | .bwt, .sa, .pac, .ann, .amb, .tar.gz | Binary; not human readable; often compressed                   |
| HISAT2 Index           | HISAT2 aligner index                        | .ht2, .ht2.tar.gz                    | Binary; not human readable; often compressed                   |
| Samtools Index         | Samtools index for FASTA                    | .fai                                 | Human readable                                                 |
| GenBank Assembly (GCA) | Archival record owned by submitter          | .gca                                 | Binary; not human readable                                     |
| RefSeq Assembly (GCF)  | NCBI-derived copy of GenBank assembly       | .gcf                                 | Binary; not human readable                                     |
| md5                    | Checksum for file integrity                 | .md5                                 | Human readable; usually plain text output from `md5sum` itself |

Would you like me to explain or elaborate on any part of this table?

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

### 🔀 Demultiplexing

Skip this section if you have already received the demultiplexed FASTQ files.

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

You can copy and paste this command if you set the variables `$run_folder`,
`$output_folder`, and `$sample_sheet` to the appropriate values.

> [!WARNING]
> As of this document's last revision, `bcl2fastq` is no longer supported;
> use `bclconvert` instead.

### 🔍 Quality Control

Quality control is an essential step in the analysis of high-throughput sequencing data.
It allows us to assess the quality of the reads and identify any issues that may
affect downstream analysis, like adapter contamination or low-quality reads.
More exciting quality issues include GC bias, mitochondrial contamination, and
over-representation of certain sequences.

| Tool                                                                     | Description                                               | Source                                                                                                                 |
| ------------------------------------------------------------------------ | --------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)[^1] | Generates html reports containing straightforward metrics | [GitHub](https://github.com/s-andrews/FastQC)                                                                          |
| [GATK](https://gatk.broadinstitute.org/hc/en-us)[^2]                     | Analyzes high-throughput sequencing data                  | [GitHub](https://github.com/broadinstitute/gatk)                                                                       |
| [Picard Tools](https://broadinstitute.github.io/picard/)[^3]             | Manipulates high-throughput sequencing data               | Comes packaged with [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) |
| [MultiQC](https://multiqc.info/)[^4]                                     | Aggregates results from bioinformatics analyses           | [GitHub](https://github.com/ewels/MultiQC)                                                                             |

[^2]: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[^3]: McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297-1303. [PMID: 20644199](https://pubmed.ncbi.nlm.nih.gov/20644199/)

[^4]: Broad Institute. Picard Tools. http://broadinstitute.github.io/picard/

[^5]: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-8. [PMID: 27312411](https://pubmed.ncbi.nlm.nih.gov/27312411/)

To run these QC applications, you need a suitable Java Runtime Environment (JRE).
Let `micromamba` handle the installation of the JRE and the tools from bioconda:

```sh
micromamba create -n qc -c conda-forge -c bioconda fastqc gatk4 picard multiqc
```

Activate the environment and run the tools:

```sh
micromamba run -n qc fastqc -o <output_dir> <fastq_file>
```

> [!TIP]
> After aligning the reads to the reference genome, these tools can be re-ran on the
> resulting SAM/BAM files to ensure that the alignment was successful or to consolidate
> the results from paired-end sequencing.

### 🏗️ Alignment

Before we can analyze the data, we need to align the reads to a reference genome.
Before aligning the reads, we need download the reference genome and build the index files.
The most recent major releases from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets)
can be found on the [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc) page.

#### Reference Genomes

<!-- TODO: add dates -->

| species                                         | assembly | release date | accession        | ftp link                                                                                 |
| ----------------------------------------------- | -------- | ------------ | ---------------- | ---------------------------------------------------------------------------------------- |
| [human](https://www.ncbi.nlm.nih.gov/grc/human) | GRCh38   | xxxx-xx-xx   | GCA_000001405.15 | [ftp](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/) |
| [mouse](https://www.ncbi.nlm.nih.gov/grc/mouse) | GRCm39   | xxxx-xx-xx   | GCA_000001635.9  | [ftp](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/)  |

> [!TIP]
> Skip building indexes from scratch and use the pre-built indexes for `bowtie2`,
> `bwa`, and `hisat2` and `samtools` in the `seqs_for_alignment_pipelines.ucsc_ids`
> folder. (It even has the 'GTT' and 'GFF' annotation files we'll need later).

#### Aligners

| Tool                                                                 | Description                                           | Key Features                                                                                                         | Best For                                        | Source                                                     |
| -------------------------------------------------------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- | ---------------------------------------------------------- |
| [BWA](https://bio-bwa.sourceforge.net/)[^6]                          | Maps short DNA sequences to reference genome          | - Uses Burrows-Wheeler Transform (BWT) for indexing<br>- Efficient for short reads<br>- Supports paired-end reads    | Whole Genome Sequencing (WGS), Exome Sequencing | [GitHub](https://github.com/lh3/bwa)                       |
| [STAR](https://github.com/alexdobin/STAR)[^7]                        | Specialized for RNA-Seq alignment                     | - Uses seed-extension search<br>- Detects novel splice junctions<br>- Fast and accurate for long reads               | RNA-Seq, especially with long reads             | [GitHub](https://github.com/alexdobin/STAR)                |
| [HISAT2](https://daehwankimlab.github.io/hisat2/)[^8]                | Splice-aware aligner for DNA and RNA sequences        | - Uses graph-based alignment<br>- Memory-efficient<br>- Supports both DNA and RNA alignment                          | RNA-Seq, WGS, particularly for large genomes    | [GitHub](https://github.com/DaehwanKimLab/hisat2)          |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)[^9] | Efficient short read aligner                          | - Uses FM-index (similar to BWT)<br>- Supports gapped, local, and paired-end alignment<br>- Memory-efficient         | ChIP-seq, WGS                                   | [GitHub](https://github.com/BenLangmead/bowtie2)           |
| [Subread](https://subread.sourceforge.net/)[^10]                     | Seed-and-vote algorithm-based aligner                 | - Fast and accurate<br>- Supports indel detection<br>- Includes read counting functionality (featureCounts)          | RNA-Seq, DNA-Seq                                | [GitHub](https://github.com/ShiLab-Bioinformatics/subread) |
| [Subjunc](https://subread.sourceforge.net/subjunc.html)[^11]         | Exon-exon junction detector (part of Subread package) | - Detects novel exon-exon junctions<br>- Uses seed-and-vote algorithm<br>- Can be used independently or with Subread | RNA-Seq, specifically for junction detection    | [GitHub](https://github.com/ShiLab-Bioinformatics/subread) |

[^6]: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168](https://pubmed.ncbi.nlm.nih.gov/19451168/)

[^7]: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21. [PMID: 23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)

[^8]: Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12(4):357-360. [PMID: 25751142](https://pubmed.ncbi.nlm.nih.gov/25751142/)

[^9]: Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357-359. [PMID: 22388286](https://pubmed.ncbi.nlm.nih.gov/22388286/)

[^10]: Liao Y, Smyth GK, Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research. 2013;41(10):e108. [PMID: 23558742](https://pubmed.ncbi.nlm.nih.gov/23558742/)

[^11]: Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research. 2019;47(8):e47. [PMID: 30783653](https://pubmed.ncbi.nlm.nih.gov/30783653/)

#### FASTQ to BAM/CRAM

(Work in progress)

### 🔬 Assembly and Quantification

Read the docs for details on experiment-specific processing and analysis:

- [RNAseq](docs/RNA-Seq.md)
- [ChIPseq](docs/ChIP-Seq.md)

## 📑 Additional Resources

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)
- [FAQs - NCBI](https://ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/)

## Acknowledgements

Shout out to these awesome docs:

- [Learn Vimscript the Hard Way](https://learnvimscriptthehardway.stevelosh.com/)
- [tao-of-tmux](https://tao-of-tmux.readthedocs.io/)
- [mamba](https://mamba.readthedocs.io/)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/index.html)
