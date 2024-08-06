library(ggplot2)

# Function to create volcano plot with labeling significant genes
create_volcano_plot <- function(results_transcripts, treatment_name) {
  # Determine significant genes based on q-value threshold (e.g., 0.05)
  results_transcripts$sig <- ifelse(results_transcripts$qval < 0.05, 
                                    ifelse(results_transcripts$log2fc > 0, "up", "down"), "notsig")
  
  # Assign colors based on significance and direction
  results_transcripts$color <- ifelse(results_transcripts$sig == "notsig", "grey", 
                                     ifelse(results_transcripts$sig == "up", "blue", "red"))
  
  # Create volcano plot
  p <- ggplot(results_transcripts, aes(x = log2fc, y = -log10(qval), color = color)) +
    geom_point() +
    scale_color_manual(values = c("grey" = "grey", "blue" = "blue", "red" = "red")) +
    theme_minimal() +
    labs(x = "Log2 Fold Change",
         y = "-log10(Q-value)",
         title = paste("Volcano Plot:", treatment_name))
  
  # Label top 10 significant genes
  top10_genes <- head(results_transcripts[order(results_transcripts$qval), ], 10)
  p <- p + geom_text(data = top10_genes, aes(label = gene_id), nudge_y = 0.5, check_overlap = TRUE)
  
  return(p)
}

# Create volcano plots for different treatments
volcano_dmso_fin <- create_volcano_plot(results_dmso_fin$transcripts, "DMSO_FIN vs Control")
volcano_dmso_oza <- create_volcano_plot(results_dmso_oza$transcripts, "DMSO_OZA vs Control")
volcano_dmso_pon <- create_volcano_plot(results_dmso_pon$transcripts, "DMSO_PON vs Control")

# Print the plots
print(volcano_dmso_fin)
print(volcano_dmso_oza)
print(volcano_dmso_pon)

