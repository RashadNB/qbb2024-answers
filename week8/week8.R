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

#Question 1
gut
# There are 13407 genes.
# They analyzed 11788 cells.
# The dimension reduction datasets present are PCA, TSNE, and UMAP

