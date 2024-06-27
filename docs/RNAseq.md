# RNAseq

Transcript-level expression analysis of RNA-seq experiments
with HISAT, StringTie and Ballgown

> [Software Availability](https://ccb.jhu.edu/software.shtml)

## Tuxedo Suite[^1]

- [HISAT2](htts://github.com/DaehwanKimLab/hisat2):
  A fast and sensitive alignment program for mapping next-generation sequencing reads
- [StringTie](https://github.com/gpertea/stringtie):
  A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts
- [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html):
  Flexible, isoform-level differential expression analysis

and optionally...

- [gffcompare](https://github.com/gpertea/gffcompare):
  A program for comparing and evaluating the accuracy of RNA-seq transcript assemblers

For more information, read the [manual](https://daehwankimlab.github.io/hisat2/manual/).

## Documentation

Protocol:
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown](https://www.nature.com/articles/nprot.2016.095)

Software:

- [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml)
- [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml)
- [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html)
- [Gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

## Workflow

### Steps 1-2: Align reads to the reference genome and sort to BAM format

See the doc on alignment [here](./docs/GENOME.md).

### Steps 3-6: `stringtie`

See the documentation in the script `stringtie.sh` for details.

### Steps 7-9: `ballgown` setup

See `ballgown.R` for details.

> [!TIP]
> Save the bg object to a file for later use:
> `save(bg, file="bg.rda")`

### Steps 10-15: Differential expression analysis

- filter to remove low-abundance genes
- identify transcripts
- identigy genes
- add gene names
- sort by p-value
- write to file

> [!IMPORTANT]
> Add log2 fold change and p-value to the gene-level data frame.
> `log2fc <- log2(gene$mean2/gene$mean1)`

[^1] Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
