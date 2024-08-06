#!/usr/bin/env Rscript
# Install core packages with BiocManager to ensure compatibility
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))


# If a package is not installed, install it and load it
install_and_load <- function(package) {
    if (!require(package, quietly = TRUE))
	BiocManager::install(package)
    library(package, character.only = TRUE)
}

install_and_load("Seurat")
install_and_load("SingleR")

# > sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.5 LTS

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

# locale:
#  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
#  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base


# loaded via a namespace (and not attached):
#   [1] Seurat_4.2.0                Rtsne_0.16                  colorspace_2.1-0            deldir_1.0-6
#   [5] ellipsis_0.3.2              ggridges_0.5.4              XVector_0.36.0              GenomicRanges_1.48.0
#   [9] rstudioapi_0.14             spatstat.data_3.0-0         leiden_0.4.3                listenv_0.8.0
#  [13] remotes_2.4.2               ggrepel_0.9.2               fansi_1.0.3                 codetools_0.2-18
#  [17] splines_4.2.1               polyclip_1.10-4             jsonlite_1.8.3              ica_1.0-3
#  [21] cluster_2.1.4               png_0.1-8                   rgeos_0.5-9                 uwot_0.1.14
#  [25] shiny_1.7.3                 sctransform_0.3.5           spatstat.sparse_3.0-0       BiocManager_1.30.18
#  [29] compiler_4.2.1              httr_1.4.4                  basilisk_1.8.1              assertthat_0.2.1
#  [33] SeuratObject_4.1.3          Matrix_1.5-1                fastmap_1.2.0               lazyeval_0.2.2
#  [37] cli_3.4.1                   later_1.3.2                 htmltools_0.5.3             tools_4.2.1
#  [41] igraph_1.3.5                GenomeInfoDbData_1.2.8      gtable_0.3.1                glue_1.7.0
#  [45] RANN_2.6.1                  reshape2_1.4.4              dplyr_1.0.10                Rcpp_1.0.9
#  [49] Biobase_2.56.0              scattermore_0.8             vctrs_0.4.2                 nlme_3.1-160
#  [53] progressr_0.11.0            zellkonverter_1.6.5         lmtest_0.9-40               spatstat.random_2.2-0
#  [57] stringr_1.5.0               globals_0.16.2              mime_0.12                   miniUI_0.1.1.1
#  [61] lifecycle_1.0.3             irlba_2.3.5.1               goftest_1.2-3               future_1.29.0
#  [65] basilisk.utils_1.8.0        zlibbioc_1.42.0             MASS_7.3-58.1               zoo_1.8-11
#  [69] scales_1.2.1                spatstat.core_2.4-4         MatrixGenerics_1.8.1        promises_1.3.0
#  [73] spatstat.utils_3.0-1        parallel_4.2.1              SummarizedExperiment_1.26.1 RColorBrewer_1.1-3
#  [77] SingleCellExperiment_1.18.1 reticulate_1.26             pbapply_1.6-0               gridExtra_2.3
#  [81] ggplot2_3.3.6               rpart_4.1.19                stringi_1.7.8               S4Vectors_0.36.2
#  [85] filelock_1.0.2              BiocGenerics_0.44.0         GenomeInfoDb_1.32.4         bitops_1.0-7
#  [89] rlang_1.0.6                 pkgconfig_2.0.3             matrixStats_1.1.0           lattice_0.20-45
#  [93] ROCR_1.0-11                 purrr_0.3.5                 tensor_1.5                  patchwork_1.1.2
#  [97] htmlwidgets_1.5.4           cowplot_1.1.1               tidyselect_1.2.0            parallelly_1.32.1
# [101] RcppAnnoy_0.0.20            plyr_1.8.8                  magrittr_2.0.3              R6_2.5.1
# [105] IRanges_2.32.0              generics_0.1.3              DelayedArray_0.22.0         DBI_1.1.3
# [109] pillar_1.8.1                mgcv_1.8-41                 fitdistrplus_1.1-8          RCurl_1.98-1.9
# [113] survival_3.4-0              abind_1.4-5                 sp_1.5-1                    dir.expiry_1.4.0
# [117] tibble_3.1.8                future.apply_1.10.0         KernSmooth_2.23-20          utf8_1.2.2
# [121] spatstat.geom_2.4-0         plotly_4.10.1               viridis_0.6.2               grid_4.2.1
# [125] data.table_1.14.4           digest_0.6.33               xtable_1.8-4                tidyr_1.2.1
# [129] httpuv_1.6.6                stats4_4.2.1                munsell_0.5.0               viridisLite_0.4.1
# Warning messages:
# 1: multiple methods tables found for 'aperm'
# 2: replacing previous import 'BiocGenerics::aperm' by 'DelayedArray::aperm' when loading 'SummarizedExperiment'
