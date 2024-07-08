# Function to create MA plot
create_ma_plot <- function(results, treatment_name) {
  ggplot(results$transcripts, aes(x = fc, y = log2fc)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Average Expression (FC)", y = "Log2 Fold Change", title = paste("MA Plot:", treatment_name))
}

ma_dmso_fin <- create_ma_plot(results_dmso_fin, "DMSO_FIN vs Control")
ma_dmso_oza <- create_ma_plot(results_dmso_oza, "DMSO_OZA vs Control")
ma_dmso_pon <- create_ma_plot(results_dmso_pon, "DMSO_PON vs Control")

print(ma_dmso_fin)
print(ma_dmso_oza)
print(ma_dmso_pon)

