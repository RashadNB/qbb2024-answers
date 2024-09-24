#!/bin/bash


# Generate a new file and give it a header separated by tabs (\t)
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt 

for i in 0.1 0.2 0.3 0.4 0.5      # for each MAF value...
    do                                                                                     
    bedtools coverage -a genome_chr1.bed -b chr1_snps_${i}.bed > snp_chr1_coverage.txt      # Find the coverage of the SNP file corresponding to i relative to chromosome 1. Assign to "snp_chr1_coverage.txt".
    snp_chr1_count=$(awk '{s+=$4}END{print s}' snp_chr1_coverage.txt)         # In the new file we made, sum up the SNPs and save as "snp_chr1_count".
    chr1_size=$(awk '{s+=$6}END{print s}' snp_chr1_coverage.txt)        # In the same file, determine the total size of chromosome 1 and save as the variable "chr1_size".
    background=$(bc -l -e "${snp_chr1_count} / ${chr1_size}")       # Save the ratio of snp_chr1_count to chr1_size as "background".

    for j in exons introns cCREs other    #For each value in the genome feature files...
        do                                                                                  
        bedtools coverage -a ${j}_chr1.bed -b chr1_snps_${i}.bed > snp_feature_coverage.txt        # Find the coverage of the SNPs in each feature and add to a new file "snp_feature_coverage.txt".
        snp_feature_count=$(awk '{s+=$4}END{print s}' snp_feature_coverage.txt)         # Add the number of SNPs in a feature.
        feature_size=$(awk '{s+=$6}END{print s}' snp_feature_coverage.txt)          # Find the number of bases in a feature.
        snp_density=$(bc -l -e "${snp_feature_count} / ${feature_size}")        # Find the ratio of counts to bases.
        enrichment=$(bc -l -e "${snp_density} / ${background}")         # Find enrichment by dividing density by background.
        echo -e "${i}\t${j}\t${enrichment}" >> snp_counts.txt       # Append the results to the file "snp_counts.txt".
        done         # End the second loop.
    done        # End the larger loop.