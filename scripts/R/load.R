# Load necessary library
library(ballgown)

# Define the data directory
data_directory <- "~/cbmf/out/ballgown_input"

# Create a vector of sample paths
sample_names <- c("Tet2RhoaG17V_24hAGX51_1", "Tet2RhoaG17V_24hAGX51_2", "Tet2RhoaG17V_24hAGX51_3",
                  "Tet2RhoaG17V_24hDMSO_1", "Tet2RhoaG17V_24hDMSO_2", "Tet2RhoaG17V_24hDMSO_3")
sample_paths <- file.path(data_directory, sample_names)

# Create phenotypic data frame
conditions <- c("AGX51", "AGX51", "AGX51", "DMSO", "DMSO", "DMSO")
pheno_data <- data.frame(sample=sample_names, condition=conditions, stringsAsFactors=FALSE)

# Load data with phenotypic data
bg <- ballgown(samples=sample_paths, pData=pheno_data, meas="all")

# Save the Ballgown object to an R data file
save(bg, file='bg.rda')

