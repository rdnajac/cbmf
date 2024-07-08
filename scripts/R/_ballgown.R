#!/usr/bin/env Rscript

# Set your directories relative to the script location
input_directory <- "./data/in/ballgown"
output_directory <- "./data/out"

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

# Example usage
agx_bg <- init_ballgown("agx")
bg_clone <- init_ballgown("clone")
dmso_fin_bg <- init_ballgown("dmso-fin")
dmso_oza_bg <- init_ballgown("dmso-oza")
dmso_pon_bg <- init_ballgown("dmso-pon")
all4groups_bg <- init_ballgown("dmso-fin-oza-pon")

# step 10: filter low abundance genes
agx_bg_filt <- subset(agx_bg, "rowVars(texpr(agx_bg)) > 1", genomesubset = TRUE)
bg_clone_filt <- subset(bg_clone, "rowVars(texpr(bg_clone)) > 1", genomesubset = TRUE)

# Filter all4groups_bg based on its own gene variance
all4groups_bg_filt <- subset(all4groups_bg, "rowVars(texpr(all4groups_bg)) > 1", genomesubset = TRUE)
# Filter dmso_fin_bg, dmso_oza_bg, dmso_pon_bg based on genes present in all4groups_bg_filt
dmso_fin_bg_filt <- subset(dmso_fin_bg, "rownames(texpr(dmso_fin_bg)) %in% rownames(texpr(all4groups_bg_filt))", genomesubset = TRUE)
dmso_oza_bg_filt <- subset(dmso_oza_bg, "rownames(texpr(dmso_oza_bg)) %in% rownames(texpr(all4groups_bg_filt))", genomesubset = TRUE)
dmso_pon_bg_filt <- subset(dmso_pon_bg, "rownames(texpr(dmso_pon_bg)) %in% rownames(texpr(all4groups_bg_filt))", genomesubset = TRUE)

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

results_agx <- run_tuxedo(agx_bg_filt, "agx")
results_clone <- run_tuxedo(bg_clone_filt, "clone")
results_dmso_fin <- run_tuxedo(dmso_fin_bg_filt, "dmso-fin")
results_dmso_oza <- run_tuxedo(dmso_oza_bg_filt, "dmso-oza")
results_dmso_pon <- run_tuxedo(dmso_pon_bg_filt, "dmso-pon")


# step 16: Identify transcripts and genes with a q value <0.05:
# subset(results_transcripts,results_transcripts$qval<0.05)
# subset(results_genes,results_genes$qval<0.05)

# step 17+: Data visualization

fpkm_boxplot <- function(bg, name) {
  fpkm <- texpr(bg, meas = 'FPKM')
  log2fpkm <- log2(fpkm + 1)
  
  samples <- colnames(pData(bg)$ids)
  groups <- unique(pData(bg)$group)
  colors <- rainbow(length(groups))
  
  # Plot boxplot with specified colors
  png(file.path(output_directory, "plots", paste0(name, "_fpkm_boxplot.png")))
  boxplot(log2fpkm, 
	  col = colors[match(pData(bg)$group, groups)],
	  las = 2,
	  ylab = "log2(FPKM+1)",
	  main = paste("FPKM for", name)
  )
  dev.off()
}
fpkm_boxplot(agx_bg, "AGX vs DMSO")
fpkm_boxplot(bg_clone, "induced vs non-induced")
fpkm_boxplot(dmso_fin_bg, "Fingolimod vs DMSO")
fpkm_boxplot(dmso_oza_bg, "Ozanimod vs DMSO")
fpkm_boxplot(dmso_pon_bg, "Ponesimod vs DMSO")
fpkm_boxplot(agx_bg_filt, "AGX vs DMSO (filtered)")
fpkm_boxplot(bg_clone_filt, "induced vs non-induced (filtered)")
fpkm_boxplot(dmso_fin_bg_filt, "Fingolimod vs DMSO (filtered)")
fpkm_boxplot(dmso_oza_bg_filt, "Ozanimod vs DMSO (filtered)")
fpkm_boxplot(dmso_pon_bg_filt, "Ponesimod vs DMSO (filtered)")

plot_pvalue_histogram <- function(pvalues, name) {
  png(file.path(output_directory, "plots", paste0(name, "_pvalues.png")))
  hist(pvalues, main = "P-value Distribution", xlab = "p-value")
  dev.off()
}

# Function to plot q-value histogram
plot_qvalue_histogram <- function(qvalues, name) {
  png(file.path(output_directory, "plots", paste0(name, "_qvalues.png")))
  hist(qvalues, main = "Q-value Distribution", xlab = "q-value")
  dev.off()
}

# Function to plot both histograms
plot_both_histograms <- function(transcripts, name) {
  plot_pvalue_histogram(transcripts$pval, name)
  plot_qvalue_histogram(transcripts$qval, name)
}

# Example usage
plot_both_histograms(results_agx$transcripts, "AGX_vs_DMSO")
plot_both_histograms(results_clone$transcripts, "Induced_vs_Non-induced")
plot_both_histograms(results_dmso_fin$transcripts, "Fingolimod_vs_DMSO")
plot_both_histograms(results_dmso_oza$transcripts, "Ozanimod_vs_DMSO")
plot_both_histograms(results_dmso_pon$transcripts, "Ponesimod_vs_DMSO")



# Function to save GCT file from ballgown object
save_gct_from_ballgown <- function(bg, prefix) {
  # Extract gene expression data
  gene_ids <- rownames(bg$expr)
  sample_ids <- colnames(bg$expr)
  expr_data <- as.matrix(bg$expr)

  # Prepare the header for the GCT file
  header <- c("Name", "Description", sample_ids)

  # Create the data frame for the GCT format
  gct_data <- data.frame(Name = gene_ids, Description = "", expr_data)
  rownames(gct_data) <- NULL

  # Write to a GCT file
  gct_file <- file.path(output_directory, paste0(prefix, "_GSEA_input.gct"))
  write.table(gct_data, file = gct_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("GCT file saved to", gct_file, "\n")
}

save_gct_from_ballgown(agx_bg, "agx")
save_gct_from_ballgown(bg_clone, "clone")
save_gct_from_ballgown(dmso_fin_bg, "dmso-fin")
save_gct_from_ballgown(dmso_oza_bg, "dmso-oza")
save_gct_from_ballgown(dmso_pon_bg, "dmso-pon")
