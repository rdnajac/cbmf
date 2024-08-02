#!/usr/bin/env Rscript
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
packages <- c(
    "Seurat", "Rtsne", "colorspace", "deldir",
    "ellipsis", "ggridges", "XVector", "GenomicRanges",
    "rstudioapi", "spatstat.data", "leiden", "listenv",
    "remotes", "ggrepel", "fansi", "codetools",
    "splines", "polyclip", "jsonlite", "ica",
    "cluster", "png", "uwot", # removing rgeos as it has been deprecated
    "shiny", "sctransform", "spatstat.sparse", "BiocManager",
    "compiler", "httr", "basilisk", "assertthat",
    "SeuratObject", "Matrix", "fastmap", "lazyeval",
    "cli", "later", "htmltools", "tools",
    "igraph", "GenomeInfoDbData", "gtable", "glue",
    "RANN", "reshape2", "dplyr", "Rcpp",
    "Biobase", "scattermore", "vctrs", "nlme",
    "progressr", "zellkonverter", "lmtest", "spatstat.random",
    "stringr", "globals", "mime", "miniUI",
    "lifecycle", "irlba", "goftest", "future",
    "basilisk.utils", "zlibbioc", "MASS", "zoo",
    "scales", "MatrixGenerics", "promises", # remove spatstat.core
    "spatstat.utils", "parallel", "RColorBrewer", "reticulate",
    "pbapply", "gridExtra", "ggplot2", "rpart",
    "stringi", "S4Vectors", "filelock", "BiocGenerics",
    "GenomeInfoDb", "bitops", "rlang", "pkgconfig",
    "matrixStats", "lattice", "ROCR", "purrr",
    "tensor", "patchwork", "htmlwidgets", "cowplot",
    "tidyselect", "parallelly", "RcppAnnoy", "plyr",
    "magrittr", "R6", "generics", "DBI",
    "pillar", "mgcv", "fitdistrplus", "RCurl",
    "survival", "abind", "sp", "dir.expiry",
    "tibble", "future.apply", "KernSmooth", "utf8",
    "spatstat.geom", "plotly", "viridis", "grid",
    "data.table", "digest", "xtable", "tidyr",
    "httpuv", "stats4", "munsell", "viridisLite"
)
for (package in packages) {
    if (!require(package, quietly = TRUE))
	BiocManager::install(package)
}
for (package in packages) {
    library(package, character.only = TRUE)
}

# now load all the packages
# answer yes when prompted to create a python environment

# detach all packages
for (package in packages) {
    detach(package, unload = TRUE)
}
# Error in detach(package, unload = TRUE) : invalid 'name' argument

# detach all packages
lapply(packages, detach, unload = TRUE)

# TODO python:
# Loading required package: GenomeInfoDb
# Unable to set up conda environment r-reticulate
# run in terminal:
# conda init
# conda create -n r-reticulate
# conda environment r-reticulate installed
# Unable to install python modules igraph and leidenalg
# run in terminal:
# conda install -n r-reticulate -c conda-forge leidenalg python-igraph pandas umap-learn
# python modules igraph and leidenalg installed
# Would you like to create a default Python environment for the reticulate package? (Yes/no/cancel)
# Would you like to create a default Python environment for the reticulate package? (Yes/no/cancel) Yes

# > sessionInfo()
# R version 4.3.3 (2024-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 24.04 LTS

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0

# locale:
#  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
# [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

# time zone: Etc/UTC
# tzcode source: system (glibc)

# attached base packages:
#  [1] grid      parallel  tools     compiler  splines   stats4    stats     graphics 
#  [9] grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] munsell_0.5.1           httpuv_1.6.15           tidyr_1.3.1            
#   [4] xtable_1.8-4            digest_0.6.36           data.table_1.15.4      
#   [7] viridis_0.6.5           viridisLite_0.4.2       plotly_4.10.4          
#  [10] utf8_1.2.4              KernSmooth_2.23-24      future.apply_1.11.2    
#  [13] tibble_3.2.1            dir.expiry_1.10.0       RCurl_1.98-1.16        
#  [16] fitdistrplus_1.2-1      survival_3.7-0          mgcv_1.9-1             
#  [19] pillar_1.9.0            DBI_1.2.3               generics_0.1.3         
#  [22] R6_2.5.1                magrittr_2.0.3          plyr_1.8.9             
#  [25] RcppAnnoy_0.0.22        parallelly_1.38.0       tidyselect_1.2.1       
#  [28] cowplot_1.1.3           htmlwidgets_1.6.4       patchwork_1.2.0        
#  [31] purrr_1.0.2             ROCR_1.0-11             lattice_0.22-6         
#  [34] pkgconfig_2.0.3         rlang_1.1.4             bitops_1.0-8           
#  [37] filelock_1.0.3          stringi_1.8.4           rpart_4.1.23           
#  [40] gridExtra_2.3           pbapply_1.7-2           RColorBrewer_1.1-3     
#  [43] spatstat.utils_3.0-5    promises_1.3.0          MatrixGenerics_1.14.0  
#  [46] matrixStats_1.3.0       scales_1.3.0            MASS_7.3-60.0.1        
#  [49] zlibbioc_1.48.2         basilisk.utils_1.14.1   future_1.34.0          
#  [52] goftest_1.2-3           irlba_2.3.5.1           lifecycle_1.0.4        
#  [55] miniUI_0.1.1.1          mime_0.12               globals_0.16.3         
#  [58] stringr_1.5.1           spatstat.random_3.3-1   spatstat.geom_3.3-2    
#  [61] spatstat.univar_3.0-0   lmtest_0.9-40           zoo_1.8-12             
#  [64] zellkonverter_1.12.1    progressr_0.14.0        nlme_3.1-165           
#  [67] vctrs_0.6.5             scattermore_1.2         Biobase_2.62.0         
#  [70] Rcpp_1.0.13             dplyr_1.1.4             reshape2_1.4.4         
#  [73] RANN_2.6.1              glue_1.7.0              gtable_0.3.5           
#  [76] GenomeInfoDbData_1.2.11 igraph_2.0.3            htmltools_0.5.8.1      
#  [79] later_1.3.2             cli_3.6.3               lazyeval_0.2.2         
#  [82] fastmap_1.2.0           assertthat_0.2.1        basilisk_1.14.3        
#  [85] reticulate_1.38.0       httr_1.4.7              BiocManager_1.30.23    
#  [88] spatstat.sparse_3.1-0   tensor_1.5              abind_1.4-5            
#  [91] sctransform_0.4.1       shiny_1.9.1             uwot_0.2.2             
#  [94] Matrix_1.6-5            png_0.1-8               cluster_2.1.6          
#  [97] ica_1.0-3               jsonlite_1.8.8          polyclip_1.10-7        
# [100] codetools_0.2-20        fansi_1.0.6             ggrepel_0.9.5          
# [103] ggplot2_3.5.1           remotes_2.5.0           listenv_0.9.1          
# [106] leiden_0.4.3.1          spatstat.data_3.1-2     rstudioapi_0.16.0      
# [109] GenomicRanges_1.54.1    GenomeInfoDb_1.38.8     XVector_0.42.0         
# [112] IRanges_2.36.0          S4Vectors_0.40.2        BiocGenerics_0.48.1    
# [115] ggridges_0.5.6          ellipsis_0.3.2          deldir_2.0-4           
# [118] colorspace_2.1-1        Rtsne_0.17              Seurat_5.1.0           
# [121] SeuratObject_5.0.2      sp_2.1-4               

# loaded via a namespace (and not attached):
#  [1] spatstat.explore_3.3-1      S4Arrays_1.2.1             
#  [3] SparseArray_1.2.4           rprojroot_2.0.4            
#  [5] RSpectra_0.16-2             here_1.0.1                 
#  [7] withr_3.0.1                 fastDummies_1.7.3          
#  [9] DelayedArray_0.28.0         rappdirs_0.3.3             
# [11] spam_2.10-0                 RcppHNSW_0.6.0             
# [13] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [15] dotCall64_1.1-1             crayon_1.5.3 


# Old version of R

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
