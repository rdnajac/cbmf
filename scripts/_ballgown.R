#!/usr/bin/env Rscript
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#' Create ballgown object from input directory.
#'
#' @param in_dir Path to directory containing ballgown input files.
#' @return Ballgown object.
#' 
#' @details
#' This function assumes that the input directory contains the following files:
#' - phenodata.csv: a CSV file containing the phenotype data.
#' in_dir/
# ├── sample01/
# │   ├── e2t.ctab
# │   ├── e_data.ctab
# │   ├── i2t.ctab
# │   ├── i_data.ctab
# │   └── t_data.ctab
# ├── sample02/
# ├── sample03/
# ├── sample04/
# ...
# └── phenodata.csv
bg_from_input_dir <- function(input_dir) {
  return(ballgown(samples = list.dirs(input_dir, full.names = TRUE, recursive = FALSE),
                  pData = read.csv(file.path(input_dir, "phenodata.csv")),
                  meas = 'all',
    )
  )
}

# Define projects and corresponding directories
data_dir <- "~/cbmf/data"

# Create ballgown objects for each project
agx_bg <- bg_from_input_dir(file.path(data_dir, "agx", "ballgown_input"))
clone_bg <- bg_from_input_dir(file.path(data_dir, "clone", "ballgown_input"))
ra_bg <- bg_from_input_dir(file.path(data_dir, "ra", "ballgown_input"))

# Save ballgown objects in current directory
save(agx_bg, file = "agx_bg.RData")
save(clone_bg, file = "clone_bg.RData")
save(ra_bg, file = "ra_bg.RData")

# Filter low abundance genes
# TODO: filter function
agx_bg_filt <- subset(agx_bg, "rowVars(texpr(agx_bg)) > 1", genomesubset = TRUE)
clone_bg_filt <- subset(clone_bg, "rowVars(texpr(clone_bg)) > 1", genomesubset = TRUE)
ra_bg_filt <- subset(ra_bg, "rowVars(texpr(ra_bg)) > 1", genomesubset = TRUE)

## DE
differential_expression_by_transcript <- function(bg) {
  return(stattest(bg, feature = 'transcript', covariate = 'group', getFC = TRUE, meas = 'FPKM')) 
}

differential_expression_by_gene <- function(bg) {
  return(stattest(bg, feature = 'gene', covariate = 'group', getFC = TRUE, meas = 'FPKM')) 
}

# function to do all that given a bg and a prefix
de_and_names <- function(bg, prefix) {
	results_transcripts <- differential_expression_by_transcript(bg) 
        results_genes <- differential_expression_by_gene(bg)
	# sanitize gene IDs
	results_genes$id <- gsub("gene-", "", results_genes$id)
	results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg), geneIDs=ballgown::geneIDs(bg), results_transcripts)
	write.csv(results_transcripts, file = paste0(prefix, "_results_transcripts.csv"), row.names = FALSE)
	write.csv(results_genes, file = paste0(prefix, "_results_genes.csv"), row.names = FALSE)
}


## Filter for genes with q-val <0.05
# subset(results_transcripts, results_transcripts$qval <=0.05)
# subset(results_genes, results_genes$qval <=0.05)

## Plotting setup
#tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

## Plotting gene abundance distribution
#fpkm <- texpr(bg_chrX, meas='FPKM')
#fpkm <- log2(fpkm +1)
#boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')

## Plot individual transcripts
#ballgown::transcriptNames(bg_chrX)[12]
#plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(bg_chrX)[12], ' : ',ballgown::transcriptNames(bg_chrX)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))

## Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

## Plot average expression
#plotMeans(ballgown::geneIDs(bg_chrX)[203], bg_chrX_filt, groupvar="sex", legend=FALSE)

