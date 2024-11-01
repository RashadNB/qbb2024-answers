library(tidyverse)
library(broom)
library(DESeq2)

#1.1
data <- read_delim("gtex_whole_blood_counts_downsample.txt")
metadata <- read_delim("gtex_metadata_downsample.txt")

# Set gene names as row names in data
data <- column_to_rownames(data, var = "GENE_NAME")

# Set subject IDs as row names in metadata
metadata <- column_to_rownames(metadata, var = "SUBJECT_ID")

#1.2
# Ensure row names in metadata match column names in count data
table(colnames(data) == rownames(metadata))

# Create DESeq2 object for differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ SEX + DTHHRDY + AGE)

#1.3
# Perform variance stabilizing transformation
vsd <- vst(dds)

# Generate and save PCA plots with grouping by SEX, DTHHRDY, and AGE
pca_plot <- function(vsd, group_var, filename) {
  p <- plotPCA(vsd, intgroup = group_var) + labs(title = paste("PCA grouped by", group_var))
  ggsave(filename = filename, plot = p, width = 6, height = 4)
}

pca_plot(vsd, "SEX", "pca_sex.png")
pca_plot(vsd, "DTHHRDY", "pca_dthhrdy.png")
pca_plot(vsd, "AGE", "pca_age.png")

# 1.3.3
# PC1 explains 48%, primarily representing cause of death; PC2 explains 7%, representing sex effects on gene expression.

#2.1
# Create data frame from vst-transformed counts
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()

# Add sample metadata to transformed data frame
vsd_df <- bind_cols(metadata, vsd_df)

# Run linear regression models for selected genes (WASH7P, SLC25A47)
run_lm <- function(gene) {
  lm(as.formula(paste(gene, "~ DTHHRDY + AGE + SEX")), data = vsd_df) %>%
    summary() %>%
    tidy()
}

m1 <- run_lm("WASH7P")
m2 <- run_lm("SLC25A47")

# WASH7P shows no significant sex-based differential expression (p = 2.79e-1).
# SLC25A47 displays significant sex-differential expression (p = 2.57e-2); males show higher expression.

#2.2
# Differential expression analysis
dds <- DESeq(dds)

# Extract sex-based differential expression results
sex_results <- results(dds, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")

# Filter and sort significant genes (FDR < 10%)
sex_results_filt <- sex_results %>%
  filter(!is.na(pvalue)) %>%
  filter(pvalue < 0.1) %>%
  arrange(pvalue)

# Count genes meeting filter criteria
sex_results_filt %>% nrow()

# 2.3.2
# 262 genes display significant differential expression between sexes.

# Add gene locations to differential expression results
fb <- read_delim("gene_locations.txt")
sex_results <- left_join(sex_results, fb, by = "GENE_NAME") %>%
  filter(!is.na(pvalue)) %>%
  arrange(pvalue)

# 2.3.3
# 17 of the top 20 genes are on the Y chromosome, with X chromosome genes also showing expected patterns.

# Retrieve differential expression data for WASH7P and SLC25A47
sex_results %>% filter(GENE_NAME == "WASH7P" | GENE_NAME == "SLC25A47")

# 2.3.4
# Results align with previous analysis: WASH7P shows no sex-differential expression; SLC25A47 does.

#2.4
# Extract differential expression results based on cause of death
death_results <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes") %>%
  as_tibble(rownames = "GENE_NAME")

# Filter and sort significant genes (FDR < 10%)
death_results <- death_results %>%
  filter(!is.na(pvalue)) %>%
  filter(pvalue < 0.1) %>%
  arrange(pvalue)

# Count genes meeting filter criteria
death_results %>% nrow()

# 2.4.1
# 16,069 genes show significant expression changes by cause of death, aligning with earlier PCA analysis.

# Create volcano plot for sex-based differential expression
ggplot(data = sex_results, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = (abs(log2FoldChange) > 1 & pvalue < 1e-1))) +
  xlim(-7, 10) +
  ylim(0, 300) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray") +
  geom_vline(xintercept = 1, linetype = 2, color = "gray") +
  geom_vline(xintercept = -1, linetype = 2, color = "gray") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "orange")) +
  labs(title = "Differential gene expression across sex",
       y = expression(-log[10]("q-value")), 
       x = expression(log[2]("fold change")))

# Save volcano plot
ggsave(filename = "volcano.png")

