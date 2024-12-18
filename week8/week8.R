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

# Question 1.1
gut
# There are 13407 genes.
# They analyzed 11788 cells.
# The dimension reduction datasets present are PCA, TSNE, and UMAP

# Question 1.2
colData(gut)
# There are 39 columns
colnames(colData(gut))

# I find the most interesting columns to me would be age, tissue, and log_n_counts.
# Together they could give some interesting insight into the fly aging process.

set.seed(88)
plotReducedDim(gut,"X_umap",color="broad_annotation")

# Question 2
# 2.1
genecounts = rowSums(assay(gut))
summary(genecounts)
# Median = 254
# Mean = 3185
# Most genes in the dataset are not highly expressed, but a few are extremely highly expressed.
sort(genecounts,TRUE)[1:3]
# Hsromega, CR45845, and roX1 are expressed the most. They are all non-coding RNAs.

# 2.2
cellcounts = colSums(assay(gut))
hist(cellcounts)
summary(cellcounts)
# Mean = 3622.
# They probably come from more common cell types.

# 2.3
celldetected = colSums(assay(gut)>0)
hist(celldetected)
summary(celldetected)
#Mean = 1059.
1059/13407 # mean over total
# This represents 7.90% of genes.

# 2.4
mito = grep("^mt:",rownames(gut),value=TRUE)
df = perCellQCMetrics(gut,subsets=list(Mito=mito))
df = as.data.frame(df)
summary(df)
colData(gut) = cbind( colData(gut), df )
plotColData(gut, y = "subsets_Mito_percent", x = "broad_annotation") +
  theme(axis.text.x = element_text(angle = 90))
#Secretory cells appear to have the most mitochondrial reads, which makes sense giving the energy-intensive role of producing and secreting proteins.

# 3.1
coi = colData(gut)$broad_annotation == "epithelial cell"
epi = gut[,coi]
plotReducedDim( epi, "X_umap", colour_by="annotation" )
marker.info = scoreMarkers( epi, colData(epi)$annotation )
chosen = marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered = chosen[order(chosen$mean.AUC, decreasing=TRUE),]
ordered$mean.AUC[1:6] #modified to show top 6 instead of top 4 like OG code.

# 3.2
# The top 6 genes are Mal-A6, Men-b, vnd, betaTry, Mal-A1, and Nhe2.
# It appears that this section of the gut is specialized in digesting carbohydrates.

plotExpression(gut, "Mal-A6", x="annotation" ) +
  theme(axis.text.x = element_text(angle = 90))

# 3.3
COI = colData(gut)$broad_annotation == "somatic precursor cell"
spc = gut[,COI]
spc.marker = scoreMarkers( spc, colData(spc)$annotation )
chosen = spc.marker[["intestinal stem cell"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
goi = rownames(ordered)[1:6]
plotExpression(spc, features = goi, x = "annotation") +
  theme(axis.text.x=element_text(angle=90)) 

# Enteroblasts and intestinal stem cells are the most similar.
# DI is the most specific marker for intestinal stem cells. 
