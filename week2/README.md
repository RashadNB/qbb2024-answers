# Part one commands

# Unzipping the downloaded file:
'tar xzf chr1_snps.tar.gz' 

# Sorting, merging and saving each of the files as a new file:
'bedtools sort -i exons.bed | bedtools merge > exons_chr1.bed'
'bedtools sort -i genes.bed | bedtools merge > genes_chr1.bed'
'bedtools sort -i cCREs.bed | bedtools merge > cCREs_chr1.bed'

# Subtracting the exons file from the genes file to get introns and saving it as a new file:
'bedtools subtract -a genes_chr1.bed -b exons_chr1.bed > introns_chr1.bed'

# Subtracting exons, introns, and regulatory sequences from the genome file to obtain everything that doesn't fall into those categories and saving it to a new file:
'bedtools subtract -a genome_chr1.bed -b exons_chr1.bed -b cCREs_chr1.bed -b introns_chr1.bed > other_chr1.bed'

# Step 2.3
- Exons have the fewest number of SNPs by far which makes sense considering the fact that they code for proteins which are more sensitive to changes that could result in deleterious effects. Rephrased, this means that they undergo the most purifying selection which can also be seen as SNP enrichment drops as MAF increases.
- A high degree of variation beyond a certain threshold is likely selected against in most populations. This process is likely being visualized in our data as a decrease in enrichment as MAFs increase since natural selection keeps below a threshold of MAF frequency. This is supported by the fact that the only genome feature that doesn't show this trend is the "other" category, which makes sense since mutations in those areas are less likely to arise as a phenotype in a population and be selected against.
- Yes. As I mentioned previously it makes sense that the coding regions would be the most conserved sequences since SNPs in those areas have the highest risk of doing damage. The non-coding areas are all grouped together which also makes sense as one would expect them to undergo similar levels of selective pressure. I initially expected the "other" feature to be higher than the rest but if SNPs in those areas impact crucial cellular components like lncRNAs it might explain why it is as conserved across the population as the others.

