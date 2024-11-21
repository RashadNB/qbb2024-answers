library(dplyr)
library(ggplot2)
library(tidyverse)

Nuclei_signals <- read.csv("Nuclei.csv")

ggplot(data = Nuclei_signals, mapping = aes(x = gene, y = nascentRNA)) + 
  geom_violin(fill = "orange") +
  labs(title = "Nascent RNA signal", 
       y = "Fluorescence",
       x = "Gene knock-down")

ggplot(data = Nuclei_signals, mapping = aes(x = gene, y = PCNA)) + 
  geom_violin(fill = "orange") +
  labs(title = "PCNA signal by gene knockdown", 
       y = "Fluorescence",
       x = "Gene knock-down")

ggplot(data = Nuclei_signals, mapping = aes(x = gene, y = log2_ratio)) + 
  geom_violin(fill = "orange") +
  labs(title = "Nascent RNA to PCNA ratio", 
       y = "Log2(nascentRNA/PCNA)",
       x = "Gene knock-down")
