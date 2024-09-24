library(tidyverse)
library(ggthemes)

snp_enrichment <- read_tsv("~/qbb2024-answers/week2/snp_counts.txt") 

snp_enrichment <- snp_enrichment %>%
  dplyr::mutate(log2_enrichment = log2(Enrichment + 1.0)) #log2 transforming enrichment.

ggplot(data = snp_enrichment, mapping = aes(x = MAF, y = log2_enrichment)) + 
  geom_line(data = snp_enrichment,
            mapping = aes(color = Feature)) + 
  labs(
    title = "SNP enrichment of genomic features",
    x = "Minor Allele Frequency",
    y = "SNP Enrichment (log2)")

ggsave(filename = "~/qbb2024-answers/week2/snp_enrichments.pdf")
