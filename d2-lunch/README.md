# day 2 lunch answers

## Answer 1
`cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c`
- This line isolates the seventh column in the document and sorts it alphabetically. It then compresses all the identical values and displays the number of each value it compressed next to each line.
`cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c | wc`
- This line counts the number of lines which reveals the number of biotypes without having to manually count all of the rows. This shows that there are 39 biotypes. The previous line shows that there are 19618 protein coding genes.
- I'd most like to learn about pseudogenes since I don't think I've ever even heard the term before and yet they make up a large number of the different biotypes.

## Answer 2
`cut -f1 hg38-gene-metadata-go.tsv | sort | uniq -c | sort -n`
- This line isolates the first column and sorts it alphabetically. Then is finds the number of times each variable in the column is represented and sorts them again numerically to reveal that the highest go_id is ENSG00000168036 with 273 entries.
`grep -w "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -k 3 -f > ENSG00000168036` 
- Here I isolate entries with the gene id of interest and sort them alphabetically based on the name_1006 column.
- Based on all the roles it seems to have in the body, this gene may code for a pretty ubiquitously expressed receptor associated with the cell membrane since it is closely associated with the Wnt pathway.

## Answer 3
`grep "IG_" genes.gtf | grep -v "_pseudogene" | cut -f 1 | uniq -c`
- This line searches for all the parts of the genes.gtf document the contain IG_ and then filters out all the parts that contain _pseudogene, leaving us with just the IG genes. It then compresses all the IG genes that are on the same chromosome and displays the number of repeats next to the chromosome name.
- 52 chr2
- 91 chr14
- 16 chr15
- 6 chr16
- 1 chr21
- 48 chr22
`grep "IG_" genes.gtf | grep "_pseudogene" | cut -f 1 | uniq -c`
- By deleting the "-v" we can isolate just the pseudogenes, showing that there are more of them and they are distributed more evenly throughout the genome.
- 1 chr1
- 45 chr2
- 1 chr8
- 5 chr9
- 1 chr10
- 84 chr14
- 6 chr15
- 8 chr16
- 1 chr18
- 48 chr22

## Answer 4
- Just using the "grep" command the way it's displayed in the question would return any line where pseudogene occurs, regardless of whether it is a "processed_pseudogene" or an "unprocessed_pseudogene" or anything esle containing the word. To fix the issue, one could employ a search algorithm using multiple filters such as the one used in Question 3.
`grep "_pseudogene" genes.gtf | grep "_processed" | cut -f9 | uniq -c`

## Answer 5
`cut -f1,4,5,14  gene-tabs.gtf > gene-tabs.bed`