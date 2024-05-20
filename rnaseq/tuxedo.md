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
if (!requireNamespace("BiocManager", quietly=TaRUE))
    install.packages("BiocManager")
BiocManager::install("ballgown")
```


## HISAT

### Genome Indexes on AWS

`aws s3 ls --no-sign-request s3://genome-idx/`

