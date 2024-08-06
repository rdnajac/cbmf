library(ggplot2)

# Function to show the distribution of FPKM values for samples colored by group based on phenodata
fpkm_boxplot <- function(bg_unfilt, bg_filt, name, output_directory = "./output") {
  fpkm_unfilt <- texpr(bg_unfilt, meas = 'FPKM')
  fpkm_filt <- texpr(bg_filt, meas = 'FPKM')
  log2fpkm_unfilt <- log2(fpkm_unfilt + 1)
  log2fpkm_filt <- log2(fpkm_filt + 1)
  
  samples_unfilt <- colnames(pData(bg_unfilt)$ids)
  samples_filt <- colnames(pData(bg_filt)$ids)
  groups <- unique(pData(bg_unfilt)$group)
  colors <- rainbow(length(groups))
  
  # Combine data into a data frame
  df_unfilt <- data.frame(log2fpkm = log2fpkm_unfilt, group = factor(pData(bg_unfilt)$group, levels = groups), condition = "Unfiltered")
  df_filt <- data.frame(log2fpkm = log2fpkm_filt, group = factor(pData(bg_filt)$group, levels = groups), condition = "Filtered")
  df_combined <- rbind(df_unfilt, df_filt)
  
  # Create ggplot object
  p <- ggplot(df_combined, aes(x = group, y = log2fpkm, fill = condition)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    labs(x = "Group", y = "log2(FPKM+1)", title = paste("FPKM for", name)) +
    theme_minimal()
  
  # Ensure the output directory exists
  if (!dir.exists(file.path(output_directory, "plots"))) {
    dir.create(file.path(output_directory, "plots"))
  }
  
  # Save the ggplot object as a PNG file using ggsave
  ggsave(file.path(output_directory, "plots", paste0(name, "_fpkm_boxplot.png")), plot = p)
}

# Example usage
fpkm_boxplot(agx_bg, agx_bg_filt, "AGX_vs_DMSO")
fpkm_boxplot(bg_clone, bg_clone_filt, "Induced_vs_Non-induced")
fpkm_boxplot(dmso_fin, dmso_fin_filt, "Fingolimod_vs_DMSO")
fpkm_boxplot(dmso_oza, dmso_oza_filt, "Ozanimod_vs_DMSO")
fpkm_boxplot(dmso_pon, dmso_pon_filt, "Ponesimod_vs_DMSO")
 
