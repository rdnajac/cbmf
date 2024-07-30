# Load necessary library
library(DESeq2)

# Step 1: Read and process count data
fc <- read.table("ra_counts.tsv", header=TRUE, row.names=1)

# Clean up column names
colnames(fc) <- gsub("_sorted_markdup.bam", "", colnames(fc))
colnames(fc) <- gsub("^\\.\\.", "", colnames(fc))

# Trim leading columns if needed (adjust indices as necessary)
fc <- fc[, 6:ncol(fc)]

# Set count data
countdata <- fc

# Step 2: Create colData DataFrame
# Define treatment vector based on your sample groups
treatment <- factor(c(
  "DMSO", "DMSO", "DMSO", 
  "Fingolimod", "Fingolimod", "Fingolimod",
  "Ozanimod", "Ozanimod", "Ozanimod",
  "Ponesimod", "Ponesimod", "Ponesimod"
))

# Create colData DataFrame
colData <- data.frame(
  row.names = colnames(countdata),  # Sample names as row names
  treatment = treatment              # Treatment conditions
)

# Step 3: Create DESeqDataSet object
ddsMat <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = colData,
  design = ~ treatment
)

# Optional: Inspect DESeqDataSet
print(ddsMat)
library(pheatmap)

# Step 4: Prepare data for heatmap
# Retrieve results from DESeq2
res <- results(dds)

# Extract top 30 genes by adjusted p-value
top_genes <- head(order(res$padj), 30)
mat <- assay(vsd)[top_genes, ]

# Center the matrix by subtracting row means
mat <- mat - rowMeans(mat)

# Create colData DataFrame with relevant columns for annotation
df <- as.data.frame(colData(vsd)[, c("cell", "treatment")])
df$cell <- factor(df$cell)
df$treatment <- factor(df$treatment)

# Generate the heatmap
pheatmap(mat, annotation_col=df, main="Heatmap of Top 30 Most Variable Genes")
# Load necessary libraries
library(pheatmap)

# Prepare the matrix for heatmap
# Extract top 30 genes by adjusted p-value
top_genes <- head(order(res$padj), 25)
mat <- assay(vsd)[top_genes, ]

# Center the matrix by subtracting row means
mat <- mat - rowMeans(mat)

# Prepare the annotation DataFrame
annotation_df <- as.data.frame(colData(vsd)[, "treatment", drop = FALSE])
colnames(annotation_df) <- "Treatment"

# Ensure the row names of annotation_df match the column names of mat
annotation_df <- annotation_df[match(colnames(mat), rownames(annotation_df)), , drop = FALSE]

# Generate the heatmap with enhanced formatting
pheatmap(
  mat,
  annotation_col = annotation_df,                # Column annotations
  main = "Top 25 Most Variable Genes", # Main title
  color = colorRampPalette(c("blue", "white", "red"))(100), # Color scheme
  cluster_rows = TRUE,                           # Cluster rows
  cluster_cols = TRUE,                           # Cluster columns
  show_rownames = TRUE,                          # Show row names
  show_colnames = TRUE,                          # Show column names
  fontsize = 10,                                 # Font size for labels
  fontsize_row = 6,                              # Font size for row names
  fontsize_col = 10,                             # Font size for column names
  border_color = NA,                             # Remove border color
  annotation_legend = TRUE,                      # Show annotation legend
  scale = "row",                                 # Scale data by row
  gaps_col = NULL,                              # Optionally add gaps between columns
  legend = TRUE                                  # Show legend
)

