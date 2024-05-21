# Tuxedo

1. HISAT
2. StringTie
3. Ballgown

## Installation

### HISAT

```sh
git clone https://github.com/DaehwanKimLab/hisat2.git && cd hisat2 && make -j $(nproc)
```

### StringTie

```sh
git clone https://github.com/gpertea/stringtie.git && cd stringtie && make release -j $(nproc)
```

### Ballgown

Start R and run:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install    ("ballgown")
```

## Genome Indexes

### AWS

`aws s3 ls --no-sign-request s3://genome-idx/`


https://genome-idx.s3.amazonaws.com/hisat/grch38_snprep.tar.gz
downlload and extract
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_snprep.tar.gz
tar -xvf grch38_snprep.tar.gz
saved in /home/ubuntu/genomes/hisat2/grch38_snprep
# add hisat2 to path
export PATH=$PATH:/home/ubuntu/src/hisat2/
