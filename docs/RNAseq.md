# RNAseq

Transcript-level expression analysis of RNA-seq experiments
with HISAT, StringTie and Ballgown

> [Software Availability](https://ccb.jhu.edu/software.shtml)

## Tuxedo Suite[^1]

[^1] Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095

Protocol:
[Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown](https://www.nature.com/articles/nprot.2016.095)

Software:

- [`HISAT2`](https://daehwankimlab.github.io/hisat2/manual/)
- [`stringtie`](https://ccb.jhu.edu/software/stringtie/index.shtml)
- [`ballgown`](https://bioconductor.org/packages/release/bioc/html/ballgown.html)
- [`gffcompare`](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

Source code:

- [HISAT2](htts://github.com/DaehwanKimLab/hisat2):
  A fast and sensitive alignment program for mapping next-generation sequencing reads
- [StringTie](https://github.com/gpertea/stringtie):
  A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts
- [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html):
  Flexible, isoform-level differential expression analysis
- [gffcompare](https://github.com/gpertea/gffcompare):
  A program for comparing and evaluating the accuracy of RNA-seq transcript assemblers

### Workflow

#### Steps 1-2: Align reads to the reference genome and sort to BAM format

See the doc on alignment.

#### Steps 3-6: `stringtie`

See the documentation in the script `stringtie.sh` for details.

### Steps 7-9: `ballgown` setup

See `_ballgown.R` for details.

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

## `> sessionInfo()`

```r
R version 4.4.0 (2024-04-24)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] ballgown_2.36.0

loaded via a namespace (and not attached):
 [1] KEGGREST_1.44.1             SummarizedExperiment_1.34.0
 [3] rjson_0.2.21                Biobase_2.64.0
 [5] lattice_0.22-6              vctrs_0.6.5
 [7] tools_4.4.0                 bitops_1.0-7
 [9] stats4_4.4.0                curl_5.2.1
[11] parallel_4.4.0              AnnotationDbi_1.66.0
[13] RSQLite_2.3.7               blob_1.2.4
[15] Matrix_1.7-0                RColorBrewer_1.1-3
[17] S4Vectors_0.42.0            GenomeInfoDbData_1.2.12
[19] compiler_4.4.0              Rsamtools_2.20.0
[21] Biostrings_2.72.1           statmod_1.5.0
[23] codetools_0.2-20            GenomeInfoDb_1.40.1
[25] RCurl_1.98-1.14             yaml_2.3.8
[27] crayon_1.5.3                BiocParallel_1.38.0
[29] DelayedArray_0.30.1         cachem_1.1.0
[31] limma_3.60.3                abind_1.4-5
[33] nlme_3.1-165                sva_3.52.0
[35] genefilter_1.86.0           locfit_1.5-9.10
[37] restfulr_0.0.15             splines_4.4.0
[39] fastmap_1.2.0               grid_4.4.0
[41] cli_3.6.3                   SparseArray_1.4.8
[43] S4Arrays_1.4.1              XML_3.99-0.17
[45] survival_3.7-0              edgeR_4.2.0
[47] UCSC.utils_1.0.0            bit64_4.0.5
[49] XVector_0.44.0              httr_1.4.7
[51] matrixStats_1.3.0           bit_4.0.5
[53] png_0.1-8                   memoise_2.0.1
[55] GenomicRanges_1.56.1        IRanges_2.38.0
[57] BiocIO_1.14.0               rtracklayer_1.64.0
[59] mgcv_1.9-1                  rlang_1.1.4
[61] Rcpp_1.0.12                 xtable_1.8-4
[63] DBI_1.2.3                   BiocManager_1.30.23
[65] BiocGenerics_0.50.0         rstudioapi_0.16.0
[67] annotate_1.82.0             jsonlite_1.8.8
[69] R6_2.5.1                    MatrixGenerics_1.16.0
[71] GenomicAlignments_1.40.0    zlibbioc_1.50.0
```

### `.Rprofile`

See [`.Rprofile`](/.Rprofile) for the code to load the libraries.

## DE across multiple groups

create multiple ballgown obects

```sh
mkdir dmso-fin dmso-oza dmso-po

for f in ./dmso-*; do cp -rv ra/ballgown/DMSO* $llgown; done
```

```r

```
