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

## 📚 Table of Contents
- [🔭 Overview](#-overview)

## 🔭 Overview

[![code style: prettier](https://img.shields.io/badge/code_style-prettier-ff69b4.svg?style=flat-square)](https://github.com/prettier/prettier)

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

Prerequisites:

- Basic command line knowledge
- A POSIX-compliant shell (e.g. `bash`, `zsh`)
- A text editor (e.g. `vim`, `nano`, `emacs`)

## Installation

Clone the repository and its submodules:

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

## 💾 Data Acquisition

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

## 📑 Additional Resources

- [Bioinformatics Stack Exchange](https://bioinformatics.stackexchange.com/)
- [Biostars](https://www.biostars.org/)
- [Biopython](https://biopython.org/)


## Acknowledgements

Shout out to these awesome docs:

- [Learn Vimscript the Hard Way](https://learnvimscriptthehardway.stevelosh.com/)
- [tao-of-tmux](https://tao-of-tmux.readthedocs.io/)
- [mamba](https://mamba.readthedocs.io/)
- [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/index.html)


<!-- [^1]: Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences. Nature Methods, 2018 doi:10.1038/s41592-018-0046-7. -->
<!-- [^2]: A. K. Maji, L. Gorenstein and G. Lentner, "Demystifying Python Package Installation with conda-env-mod," 2020 IEEE/ACM International Workshop on HPC User Support Tools (HUST) and Workshop on Programming and Performance Visualization Tools (ProTools), GA, USA, 2020, pp. 27-37, doi: 10.1109/HUSTProtools51951.2020.00011. -->
