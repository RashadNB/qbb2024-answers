library("tidyverse")
df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 )
df %>%
  group_by(SUBJECT) %>%
  summarize(count = n()) %>% #summarize the SUBJECT data by number of subjects.
  arrange(count)       #arranges data from lowest to highest.
# GTEX-1JMI6 and GTEX-1PAR6 have the smallest number of subjects, with 1 in each.

df %>%
  group_by(SUBJECT) %>%
  summarize(count = n()) %>%
  arrange(-count) #same as the code above but the "-" reverses the order.
# K-562 has the most subjects at 217.

df %>%
  group_by(SMTSD) %>%
  summarize(count2 = n()) %>%
  arrange(count2)
# Kidney has the least samples at 4 and whole blood has the most at 3288. Tissue diversity and accessibility might play a role. The brain, for instance, performs myriad functions with diverse tissues and structures whereas the kidneys may be relatively homogenous.

df_npj8 <- df %>%
  filter(SUBJECT == "GTEX-NPJ8") #assigning a single subject's data as a new object.

df_npj8  %>%
  group_by(SMTSD) %>%
  summarize(count3 = n()) %>%
  arrange(-count3)
# Many tissues are tied at 1 sample and are listed alphabetically.
# Whole blood samples are by far the most numerous at 9.

view(df_npj8)
#The sequencing methods are different, as well as the method of isolating the nucleotides. Some use DNA isoplation and some are RNA isolation.

SMATSSCR <- df %>%
  filter( !is.na(SMATSSCR) ) %>%
  group_by(SUBJECT) %>%
  summarize(meanSM=mean(SMATSSCR))

  sum(SMATSSCR$meanSM == 0) #There are 15 subjects with a mean score of 0
  
# I would be curious to identify the subjects with a mean score of 0 and determine if there are any relevant patterns between them like tissue type.
# Using ggplot to visualize the data in the form of a boxplot would be an effective way of presenting it in a report. Generally displaying data visually allows it to be more digestible and have conclusions drawn from it more effectively.