# Load the necessary libraries
library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(ggfortify)
library(ggplot2)

# Check session information
sessionInfo()

# Load the data using read_tsv
data <- readr::read_tsv("salmon.merged.gene_counts.tsv")

# Set the gene_name column as row names
data <- column_to_rownames(data, var = "gene_name")

# Remove the 'gene_id' column
data <- data %>% dplyr::select(-gene_id)

# Convert numeric columns to integers using mutate_if
data <- data %>% dplyr::mutate_if(is.numeric, as.integer)

# Filter rows to keep only genes with at least 100 reads across all samples
data <- data[rowSums(data) > 100, ]

# Select narrow region samples (midgut sections)
narrow_data <- dplyr::select(data, "A1_Rep1":"P2-4_Rep3")

# Generate tibble with tissues and replicate numbers based on sample names
narrow_metadata <- tibble(Tissue = as.factor(c("A1", "A1", "A1",
                                               "A2-3", "A2-3", "A2-3",
                                               "Cu", "Cu", "Cu",
                                               "LFC-Fe", "LFC-Fe", "Fe",
                                               "LFC-Fe", "Fe", "Fe",
                                               "P1", "P1", "P1",
                                               "P2-4", "P2-4", "P2-4")),
                         Replicate = as.factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 1, 3, 2, 3, 1, 2, 3, 1, 2, 3)))

# Create DESeqDataSet object
narrow <- DESeqDataSetFromMatrix(countData = as.matrix(narrow_data), 
                                 colData = narrow_metadata, design = ~Tissue)

# Run DESeq function to fit the model
dds_narrow <- DESeq(narrow)

# Perform variance-stabilizing transformation (VST)
vstnarrow <- vst(dds_narrow)

# Plot mean-variance
meanSdPlot(assay(vstnarrow))

# PCA plot of the VST-transformed data
narrowPcaData <- plotPCA(vstnarrow, intgroup = c("Replicate", "Tissue"), returnData = TRUE)

# Plot PCA with ggplot
ggplot(narrowPcaData, aes(PC1, PC2, color = Tissue, shape = Replicate)) +
  geom_point(size = 5) +
  labs(title = "PCA of Narrow Region Samples",
       x = paste0("PC1: ", round(100 * attr(narrowPcaData, "percentVar")[1], 2), "% variance"),
       y = paste0("PC2: ", round(100 * attr(narrowPcaData, "percentVar")[2], 2), "% variance")) +
  theme_minimal()

#convert to matrix
narrowmatrix = as.matrix(assay(vstnarrow))

# Find means one set of replicates at a time
combined = narrowmatrix[,seq(1, 21, 3)]
combined = combined + narrowmatrix[,seq(2, 21, 3)]
combined = combined + narrowmatrix[,seq(3, 21, 3)]
combined = combined / 3

# Set seed for reproducibility
set.seed(8)

# Perform k-means clustering on the VST-transformed data (12 clusters)
k_clusters <- kmeans(assay(vstnarrow), centers = 12)$cluster

# Order genes based on the cluster assignments
ordering <- order(k_clusters)

# Reorder clusters to match the gene ordering
k_clusters <- k_clusters[ordering]

# Plot heatmap of expression data with genes ordered by cluster
heatmap(assay(vstnarrow)[ordering, ], 
        Rowv = NA, 
        Colv = NA, 
        RowSideColors = RColorBrewer::brewer.pal(12, "Paired")[k_clusters])

# Extract gene names from cluster 1
cluster_1_genes <- rownames(assay(vstnarrow)[k_clusters == 1, ])

# Save the gene names from cluster 1 to a text file (for use with PANTHER)
write.table(cluster_1_genes, "cluster_genes.txt", sep = "\n", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
