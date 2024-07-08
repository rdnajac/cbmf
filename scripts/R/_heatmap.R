# Function to create heatmap
create_heatmap <- function(results, plot_title) {
  # Assuming you want to visualize top 100 DE genes based on q-value
  top_de_genes <- results$transcripts[order(results$transcripts$qval), ][1:100, ]
  
  mat <- as.matrix(top_de_genes$log2fc)
  rownames(mat) <- top_de_genes$gene_name
  
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE,
           annotation_col = top_de_genes$qval,
           main = plot_title)
}

heatmap_dmso_fin <- create_heatmap(results_dmso_fin, "Heatmap of Top DE Genes: DMSO_FIN vs Control")
heatmap_dmso_oza <- create_heatmap(results_dmso_oza, "Heatmap of Top DE Genes: DMSO_OZA vs Control")
heatmap_dmso_pon <- create_heatmap(results_dmso_pon, "Heatmap of Top DE Genes: DMSO_PON vs Control")

print(heatmap_dmso_fin)
print(heatmap_dmso_oza)
print(heatmap_dmso_pon)

