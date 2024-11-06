library(tidyverse)
library(DropletUtils)
library(zellkonverter)
library(ggthemes)
library(scuttle)
library(scater)
library(scran)

gut = readH5AD("v2_fca_biohub_gut_10x_raw.h5ad")
assayNames(gut) = "counts"
gut = logNormCounts(gut)

#Question 1.1
gut
# There are 13407 genes.
# They analyzed 11788 cells.
# The dimension reduction datasets present are PCA, TSNE, and UMAP

#Question 1.2
colData(gut)
# There are 39 columns
colnames(colData(gut))

# I find the most interesting columns to me would be age, tissue, and log_n_counts.
# Together they could give some interesting insight into the fly aging process.

set.seed(88)
plotReducedDim(gut,"X_umap",color="broad_annotation")

