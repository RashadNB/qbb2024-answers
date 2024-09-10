library(tidyverse)

expression <- read_tsv("qbb2024-answers/d4/dicts_expr.tsv")

glimpse(expression)

expression <- expression %>%
  mutate(Tissue_Data = paste0(Tissue, " ", GeneID)) %>%
  mutate(Log2_Expr = log2(Expr + 1))

ggplot(data = expression,
       mapping = aes(x = Tissue_Data, y = Log2_Expr)) +
  geom_violin() +
  coord_flip() +
  labs(x = "Tissue Type + Gene", y = "Log2 Gene Expression")
