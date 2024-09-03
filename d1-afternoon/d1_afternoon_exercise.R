library(tidyverse)
df <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
glimpse(df)
df_RNA <- df %>% filter(SMGEBTCHT == "TruSeq.v1")
view(df_RNA)
# Questions 1-3


ggplot(data = df_RNA,
       mapping = aes( x = SMTSD)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# Question 4

ggplot(data = df_RNA, mapping = aes( x = SMRIN)) +
  geom_histogram(binwidth = 0.1) #Distribution appears biomodal, with most samples clustered between 6 and 8 and others being close to 10
# Question 5

ggplot(data = df_RNA, mapping = aes( x = SMTSD, y = SMRIN)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Question 6  All the highest quality samples come from cultured cells. It is probably easier to extract RNA, or redo extractions, with a large cultured cell population than with pdx samples.
# Other than that the only obvious differences between tissues are with the sample size. Notably, kidney medulla samples are scarce and lower quality.

ggplot(data = df_RNA, mapping = aes( x = SMTSD, y = SMGNSDTC)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Question 7.  The testis and whole blood stand out as outliers. The testis seem to express more genes and whole blood express fewer genes.
# The provided paper posits that male germline cells express many genes in order to correct mutations and maintain the integrity of the genome.
# Whole blood may express fewer genes since those cells are not actively replicating for the most part, although this is not covered by the paper.

ggplot(data = df_RNA, mapping = aes( x = SMTSISCH, y = SMRIN)) +
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(. ~ SMTSD) +
  geom_smooth(method = lm)
#Question 8
# For most tissues RNA quality degrades as ischemic time increases. This effect is less pronounced for tissues like most of the areas of the brain.

ggplot(data = df_RNA, mapping = aes( x = SMTSISCH, y = SMRIN)) +
  geom_point(size = 0.5, alpha = 0.5, mapping = aes(color = SMATSSCR)) +
  facet_wrap(. ~ SMTSD) +
  geom_smooth(method = lm)
#Question 9
# Autolysis score increases with ischemic time, which in turn results in lower quality RNA.
# Some tissues display no autoloysis score, indicating that perhaps they do not autolyse under the conditions they were subjected to here. In those lines RNA quality stayed relatively constant.

ggplot(data = df_RNA, mapping = aes( x = SMNABTCHT, y = SMRIN)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Question 10 part 1
# Trizol extraction seems to preserve RNA quality the best.
ggplot(data = df_RNA, mapping = aes( x = SMNABTCHT, y = SMRIN)) +
  geom_boxplot() +
  facet_wrap(. ~ SMTSD) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Question 10 part 2
# Plotting again by cell type shows how part 1 alone could be misleading. Only the cultured cancer cells can be extracted using Trizol which may be why that method appears to be superior.
