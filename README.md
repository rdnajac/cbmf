# 🧬 Combinatorial Bioinformatic Meta-Framework

A single point of access to thousands of biomedical research software packages.

## 🔭 Overview

The Combinatorial Bioinformatic Meta-Framework (CBMF) is a collection of
tools and resources for analyzing high-throughput sequencing data. It
leverages the power of the `micromamba` package manager to install software
and manage dependencies. The framework is designed to be modular, allowing
users to select the tools they need for their specific analyses.

## 🚀 Quick Start

Clone the repository...

```sh
git clone https://github.com/rdnajac/cbmf
```

and run initialization script.

```sh
cd cbmf && $SHELL <INSTALL
```

Then skip ahead to the [Sequencing Workflows](#-sequencing-workflows) section.
Otherwise, read on for more information about the CBMF, how it works,
and how to use it.

## 📚 Documentation

CBMF is designed to be user-friendly, with detailed documentation and tutorials
to help users get started. The documentation is organized into sections based
on the type of analysis being performed, such as RNA-Seq, ChIP-Seq, and WGS.

There is also a Jira-style wiki with additional resources and tutorials.
You can access the wiki by clicking on the tab at the top of the page
or by following [this link](https://github.com/rdnajac/cbmf/wiki).

> [!TIP]
> These files are written using [GitHub Flavored Markdown \(GFM\)](https://github.github.com/gfm/),
> a superset of the [original](https://daringfireball.net/projects/markdown/syntax)
> lightweight markup language with plain text formatting syntax.
> GFM gives us additional features (like tables, footnotes, and  
> like this) while maintaining readability in plain text.

## 📦 Package Management

Bundled with the CBMF is the lightweight package manager
[`micromamba`](https://mamba.readthedocs.io/en/latest/index.html)
with access to the entire suite of bioinformatics software available on
[Bioconda](https://bioconda.github.io/index.html)[^bioconda], a channel for the
[conda](https://docs.conda.io/en/latest/) package manager (including
all available Bioconductor[^bioconductor] software.

Micromamba is not a `conda` distribution, but a statically linked C++ executable
that can be used to install `conda` environments. It is a lightweight binary
that handles the installation of conda environments without root privileges,
or the need for a base environment or a Python installation, making it ideal
for use in high-performance computing clusters.

If you want to install it on your own, skip the init scripts and run:

```
"$SHELL" <(curl -L micro.mamba.pm/install.sh)
```

CBMF comes with some sensible defaults and pre-configured environments for
common bioinformatics tasks, but users can easily create their own environments
using the [Mamba API](https://mamba.readthedocs.io/en/latest/index.html).

## 🧬 Sequencing Workflows

The following sections outline the steps involved in analyzing
high-throughput sequencing data:

1. [Demultiplexing](#-demultiplexing)
2. [Quality Control](#-quality-control)
3. [Alignment](#-alignment)
4. [Assembly and Quantification](#-assembly-and-quantification)

But first, a quick note on file types and formats:

> [!CAUTION]
> This table is a work in progress and may not be complete.

| File Type | Description                                 | Typical Extension |
| --------- | ------------------------------------------- | ----------------- |
| FASTQ     | Raw reads from sequencer (often compressed) | .fastq, .fastq.gz |
| FASTA     | Sequence data                               | .fasta, .fa, .fna |
| GTF       | Gene Transfer Format (GFF2 variant)         | .gtf, .gtf.gz     |
| GFF       | General Feature Format                      | .gff3, .gff3.gz   |
| SAM       | Sequence Alignment/Map                      | .sam, .sam.gz     |
| BAM       | Binary Alignment/Map                        | .bam              |
| CRAM      | Compressed Reference-oriented Alignment Map | .cram             |
| VCF       | Variant Call Format                         | .vcf,             |
| BED       | Browser Extensible Data                     | .bed,             |
| TSV       | Tab-Separated Values                        | .tsv,             |
| CSV       | Comma-Separated Values                      | .csv,             |
| JSON      | JavaScript Object Notation                  | .json,            |
| Markdown  | Markup language for documentation           | .md, .markdown    |
| TXT       | Plain Text                                  | .txt,             |
| md5       | Checksum for file integrity                 | .md5, .txt        |

### 🔀 Demultiplexing

Skip this section and read about [Quality Control](#-quality-control)
if you have already received the demultiplexed FASTQ files.

The following command is the default [`bcl2fastq`](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
command for demultiplexing on the Nextseq, but with the `--no-lane-splitting`
option added to combine the reads from all four lanes into a single FASTQ file:

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
> use [`bclconvert`](https://emea.support.illumina.com/downloads/bcl-convert-user-guide.html)
> if you have a used recent Illumina sequencer (NovaSeq, NextSeq 1000/2000, etc.).

### 🔍 Quality Control

Quality control is an essential step in the analysis of high-throughput sequencing
data. It allows us to assess the quality of the reads and identify any issues
that may affect downstream analysis, like adapter contamination or low-quality
reads. More interesting quality issues include GC bias, mitochondrial
contamination, and over-representation of certain sequences.

| Tool                                                                          | Description                                               | Source                                                                                                                 |
| ----------------------------------------------------------------------------- | --------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)[^fastqc] | Generates html reports containing straightforward metrics | [GitHub](https://github.com/s-andrews/FastQC)                                                                          |
| [GATK](https://gatk.broadinstitute.org/hc/en-us)[^gatk]                       | Analyzes high-throughput sequencing data                  | [GitHub](https://github.com/broadinstitute/gatk)                                                                       |
| [Picard Tools](https://broadinstitute.github.io/picard/)                      | Manipulates high-throughput sequencing data               | Comes packaged with [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) |
| [MultiQC](https://multiqc.info/)[^4]                                          | Aggregates results from bioinformatics analyses           | [GitHub](https://github.com/ewels/MultiQC)                                                                             |

To run these QC applications, you need a suitable Java Runtime Environment (JRE).
Let `micromamba` handle the installation of the JRE and the tools from bioconda:

```sh
micromamba create -n qc -c conda-forge -c bioconda fastqc gatk4 picard multiqc
micromamba run -n qc fastqc -o <output_dir> <fastq_file>
```

You can also use the qc configuration file in the `dev` folder to create the environment:

```sh
micromamba env create -f dev/qc.yml
micromamba activate qc
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

| Tool                                                                       | Description                                           | Key Features                                                                                                         | Best For                                        | Source                                                     |
| -------------------------------------------------------------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- | ---------------------------------------------------------- |
| [BWA](https://bio-bwa.sourceforge.net/)[^bwa]                              | Maps short DNA sequences to reference genome          | - Uses Burrows-Wheeler Transform (BWT) for indexing<br>- Efficient for short reads<br>- Supports paired-end reads    | Whole Genome Sequencing (WGS), Exome Sequencing | [GitHub](https://github.com/lh3/bwa)                       |
| [STAR](https://github.com/alexdobin/STAR)[^star]                           | Specialized for RNA-Seq alignment                     | - Uses seed-extension search<br>- Detects novel splice junctions<br>- Fast and accurate for long reads               | RNA-Seq, especially with long reads             | [GitHub](https://github.com/alexdobin/STAR)                |
| [HISAT2](https://daehwankimlab.github.io/hisat2/)[^hisat2]                 | Splice-aware aligner for DNA and RNA sequences        | - Uses graph-based alignment<br>- Memory-efficient<br>- Supports both DNA and RNA alignment                          | RNA-Seq, WGS, particularly for large genomes    | [GitHub](https://github.com/DaehwanKimLab/hisat2)          |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)[^bowtie2] | Efficient short read aligner                          | - Uses FM-index (similar to BWT)<br>- Supports gapped, local, and paired-end alignment<br>- Memory-efficient         | ChIP-seq, WGS                                   | [GitHub](https://github.com/BenLangmead/bowtie2)           |
| [Subread](https://subread.sourceforge.net/)[^subread]                      | Seed-and-vote algorithm-based aligner                 | - Fast and accurate<br>- Supports indel detection<br>- Includes read counting functionality (featureCounts)          | RNA-Seq, DNA-Seq                                | [GitHub](https://github.com/ShiLab-Bioinformatics/subread) |
| [Subjunc](https://subread.sourceforge.net/subjunc.html)[^subjunc]          | Exon-exon junction detector (part of Subread package) | - Detects novel exon-exon junctions<br>- Uses seed-and-vote algorithm<br>- Can be used independently or with Subread | RNA-Seq, specifically for junction detection    | [GitHub](https://github.com/ShiLab-Bioinformatics/subread) |

#### FASTQ to BAM/CRAM

(Work in progress)

### 🔬 Assembly and Quantification

Read the wiki for details on experiment-specific processing and analysis.

## 📑 Additional Resources

Writing guides:

- [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).
- [About READMEs](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes)
- [Art of README](https://github.com/hackergrrl/art-of-readme)

Nerd stuff:

- [Write a Good Technical Report](https://ieeexplore.ieee.org/document/6448763)
- [Code Documentation](https://ieeexplore.ieee.org/abstract/document/5484109)
- [Semantic line breaks](https://sembr.org/)
- [Semantic Versioning](https://semver.org/)

Message boards:

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)

FAQs:

- [FAQs - NCBI](https://ncbi.nlm.nih.gov/datasets/docs/v2/troubleshooting/faq/)

## 👍 Acknowledgements

Shout out to these awesome docs:

- [Learn Vimscript the Hard Way](https://learnvimscriptthehardway.stevelosh.com/)
- [tao-of-tmux](https://tao-of-tmux.readthedocs.io/)
- [mamba](https://mamba.readthedocs.io/)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/index.html)

Thank you to my labmates in the [Palomero Lab](http://palomerolab.org/)
for their feedback and guidance.

<!-- Refereneces -->

[^bioconda]:
    Grüning, B., Dale, R., Sjödin, A. et al. Bioconda: sustainable and comprehensive
    software distribution for the life sciences. Nat Methods 15, 475–476 (2018).
    https://doi.org/10.1038/s41592-018-0046-7

[^bioconductor]:
    Gentleman, R.C., Carey, V.J., Bates, D.M. et al. Bioconductor: open software
    development for computational biology and bioinformatics.
    Genome Biol 5, R80 (2004). https://doi.org/10.1186/gb-2004-5-10-r80

[^fastqc]:
    Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
    Available online at: <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

[^gatk]:
    McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce
    framework for analyzing next-generation DNA sequencing data.
    Genome Res. 2010;20(9):1297-1303. [PMID: 20644199](https://pubmed.ncbi.nlm.nih.gov/20644199/)

[^bwa]:
    Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform.
    Bioinformatics, 25:1754-60. [PMID: 19451168](https://pubmed.ncbi.nlm.nih.gov/19451168/)

[^star]:
    Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner.
    Bioinformatics. 2013;29(1):15-21. [PMID: 23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)

[^hisat2]:
    Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements.
    Nat Methods. 2015;12(4):357-360. [PMID: 25751142](https://pubmed.ncbi.nlm.nih.gov/25751142/)

[^bowtie2]:
    Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2.
    Nat Methods. 2012;9(4):357-359. [PMID: 22388286](https://pubmed.ncbi.nlm.nih.gov/22388286/)

[^subread]:
    Liao Y, Smyth GK, Shi W. The Subread aligner: fast, accurate and scalable
    read mapping by seed-and-vote. Nucleic Acids Research. 2013;41(10):e108. [PMID: 23558742](https://pubmed.ncbi.nlm.nih.gov/23558742/)

[^subjunc]:
    Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and
    better for alignment and quantification of RNA sequencing reads.
    Nucleic Acids Research. 2019;47(8):e47. [PMID: 30783653](https://pubmed.ncbi.nlm.nih.gov/30783653/)
