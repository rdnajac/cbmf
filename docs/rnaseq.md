# RNAseq

## Tuxedo Suite\[^1\]

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
> This function makes assumptions about the software being installed.
> It may not work for all software.

```bash
install_and_add_to_path https://github.com/DaehwanKimLab/hisat2.git
install_and_add_to_path https://github.com/gpertea/stringtie.git
```

### HISAT2

To install, run this in the folder containing the executables:

```bash
sudo cp -v hisat2 hisat2-align-s hisat2-align-l hisat2-build hisat2-build-s hisat2-build-l hisat2-inspect hisat2-inspect-s hisat2-inspect-l /usr/local/bin
```


```bash

### Ballgown

Start R and run:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ballgown")
```


