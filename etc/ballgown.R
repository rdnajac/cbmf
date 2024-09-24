#!/usr/bin/env Rscript

# Set your directories relative to the script location
input_directory <- "./"
output_directory <- "./"

#' Initialize or load ballgown object from input directory.
#'
#' @param project_directory Sub-directory within input_directory containing ballgown input files.
#' @return Ballgown object.
#' 
#' @details
#' This function assumes that the project directory contains the following files:
#' ├── sample01/
#' │   ├── e2t.ctab
#' │   ├── e_data.ctab
#' │   ├── i2t.ctab
#' │   ├── i_data.ctab
#' │   └── t_data.ctab
#' ├── sample02/
#' ├── sample03/
#' ├── sample04/
#' ...
#' └── phenodata.csv (ids,group\nsample01,group1\nsample02,group2\n...)
init_ballgown <- function(project_directory) {
  bg_rda <- file.path(output_directory, paste0(project_directory, "_bg.RData"))

  if (!file.exists(bg_rda)) {
    # Create new ballgown object if it doesn't exist
    sample_directory <- file.path(input_directory, project_directory)
    samples_list <- list.dirs(sample_directory, full.names = TRUE, recursive = FALSE)
    pdata_from_existing_csv <- read.csv(file.path(sample_directory, "phenodata.csv"))
    bg <- ballgown(samples = samples_list, pData = pdata_from_existing_csv, meas = 'all')
    cat("Saving ballgown object to", bg_rda, "\n")
    save(bg, file = bg_rda)
  } else {
    cat("Loading existing ballgown object from", bg_rda, "\n")
    load(bg_rda)
  }
  return(bg)
}


# samples_list <- list.dirs(sample_directory, full.names = TRUE, recursive = FALSE)
# again, but just for this cwd
samples_list <- list.dirs(".", full.names = TRUE, recursive = FALSE)
# ignore hidden

# Example usage
# agx_bg <- init_ballgown("agx")
# bg_clone <- init_ballgown("clone")
dmso_fin_bg <- init_ballgown("dmso-fin")
dmso_oza_bg <- init_ballgown("dmso-oza")
dmso_pon_bg <- init_ballgown("dmso-pon")
# all4groups_bg <- init_ballgown("dmso-fin-oza-pon")

# step 10: filter low abundance genes
# https://www.biostars.org/p/367562/
install.packages('metaMA')
library(metaMA)
?rowVars
# >bg_chrX_filt = subset(bg_chrX,″rowVars(texpr(bg_chrX)) >1″,genomesubset=TRUE)
dmso_fin_bg_filt <- subset(dmso_fin_bg, "rowVars(texpr(dmso_fin_bg)) > 1", genomesubset = TRUE)
dmso_oza_bg_filt <- subset(dmso_oza_bg, "rowVars(texpr(dmso_oza_bg)) > 1", genomesubset = TRUE)
dmso_pon_bg_filt <- subset(dmso_pon_bg, "rowVars(texpr(dmso_pon_bg)) > 1", genomesubset = TRUE)

# function to do all that given a bg and a prefix
run_tuxedo <- function(bg, prefix) {

	# step 11: differential expression analysis of transcripts
	results_transcripts <- stattest(bg, feature = 'transcript', 
					    covariate = 'group', 
					    getFC = TRUE, 
					    meas = 'FPKM')

	# step 12: differential expression analysis of genes
        results_genes <- stattest(bg, feature = 'gene',
				      covariate = 'group',
				      getFC = TRUE,
				      meas = 'FPKM')

	# sanitize gene IDs
	results_genes$id <- gsub("gene-", "", results_genes$id)

	# step 13: add gene names and ids to the results
	results_transcripts <- data.frame(gene_name=ballgown::geneNames(bg), 
					  gene_id=ballgown::geneIDs(bg), 
					  results_transcripts)

	# step 14: sort the results from smallest p-value to largest
	results_transcripts <- arrange(results_transcripts, pval)
	results_genes <- arrange(results_genes, pval)

	# # sort by q-value instead
	# results_transcripts <- arrange(results_transcripts, qval)
	# results_genes <- arrange(results_genes, qval)

	# Calculate log2 fold change
	results_genes$log2fc <- log2(results_genes$fc +1)
	results_transcripts$log2fc <- log2(results_transcripts$fc +1)

	# remove the 'feature' column prior to saving
	results_transcripts$feature <- NULL
	results_genes$feature <- NULL

	# step 15: save the results to a file (and return results object for further analysis)
	write.csv(results_transcripts, file = file.path(output_directory, paste0(prefix, "_DE_results_transcripts.csv")), row.names = FALSE)
	write.csv(results_genes, file = file.path(output_directory, paste0(prefix, "_DE_results_genes.csv")), row.names = FALSE)
	cat("Results saved to", file.path(output_directory, paste0(prefix, "_results_transcripts.csv")), "\n")
	return(list(transcripts = results_transcripts, genes = results_genes))
}
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
results_dmso_fin <- run_tuxedo(dmso_fin_bg_filt, "dmso-fin")
results_dmso_oza <- run_tuxedo(dmso_oza_bg_filt, "dmso-oza")
results_dmso_pon <- run_tuxedo(dmso_pon_bg_filt, "dmso-pon")


# step 16: Identify transcripts and genes with a q value <0.05:
# subset(results_transcripts,results_transcripts$qval<0.05)
# subset(results_genes,results_genes$qval<0.05)

fin_sig_transcripts <- subset(results_dmso_fin$transcripts, results_dmso_fin$transcripts$qval < 0.05)
fin_sig_genes <- subset(results_dmso_fin$genes, results_dmso_fin$genes$qval < 0.05)
oza_sig_transcripts <- subset(results_dmso_oza$transcripts, results_dmso_oza$transcripts$qval < 0.05)
oza_sig_genes <- subset(results_dmso_oza$genes, results_dmso_oza$genes$qval < 0.05)
pon_sig_transcripts <- subset(results_dmso_pon$transcripts, results_dmso_pon$transcripts$qval < 0.05)
pon_sig_genes <- subset(results_dmso_pon$genes, results_dmso_pon$genes$qval < 0.05)

# pretty print the results
# top_DE_transcripts_by_treatment_group_v_dmso <- data.frame(fin=fin_sig_transcripts$gene_name, oza=oza_sig_transcripts$gene_name, pon=pon_sig_transcripts$gene_name)

  # arguments imply differing number of rows: 3, 1, 0
# make it so that they can have empty rows
# top_DE_transcripts_by_treatment_group_v_dmso <- data.frame(fin=fin_sig_transcripts$gene_name, oza=oza_sig_transcripts$gene_name, pon=pon_sig_transcripts$gene_name)

# DE_transcripts_fin <- results_dmso_fin$transcripts
# sort by q-value in the same step
DE_transcripts_fin <- results_dmso_fin$transcripts[order(results_dmso_fin$transcripts$qval),]
DE_transcripts_oza <- results_dmso_oza$transcripts[order(results_dmso_oza$transcripts$qval),]
DE_transcripts_pon <- results_dmso_pon$transcripts[order(results_dmso_pon$transcripts$qval),]

# save to csv
write.csv(DE_transcripts_fin, file = file.path(output_directory, "DE_transcripts_fin_dmso.csv"), row.names = FALSE)
write.csv(DE_transcripts_oza, file = file.path(output_directory, "DE_transcripts_oza_dmso.csv"), row.names = FALSE)
write.csv(DE_transcripts_pon, file = file.path(output_directory, "DE_transcripts_pon_dmso.csv"), row.names = FALSE)


