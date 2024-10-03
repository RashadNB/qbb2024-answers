#!/usr/bin/env bash 

### Question 2.1 ###
grep -c "^>" sacCer3.fa # grep for > because chromosomes start with > in FASTA files.
# There are 17 chromosomes in the yeast genome (16 plus one mitochondrial).

### Question 2.2 ###
for my_sample in A01_*.fastq # align reads to reference genome by extracting the basename in each seq file, align it to the reference, and put it in a sam file.
do
    echo ${my_sample}
    my_sample=$(basename ${my_sample} .fastq)
    bwa mem -t 4 -R "@RG\tID:${my_sample}\tSM:${my_sample}" sacCer3.fa ${my_sample}.fastq > ${my_sample}.sam
done
less -S A01_09.sam
grep -v "^@" A01_09.sam | wc -l 
# There are 669548 read alignments in the file.


### Question 2.3 ###
grep -w "chrIII" A01_09.sam | wc -l # Finds lines in the sam file that reference chromosome 3.
# There are 18196 alignments to ch3.


### Question 2.4 ###
for my_sample in A01_*.sam
do
    my_sample=$(basename ${my_sample} .sam)
    samtools sort -@ 4 -O bam -o ${my_sample}.bam ${my_sample}.sam
    samtools index ${my_sample}.bam
done
# Coverage does match, but almost on a technicality because of how not uniform it is. Some areas have disproprtionately more coverage than others so an "average" read depth may not be representative.


### Question 2.5 ###
# There are 3 SNPs visible. One of them is only convered by two reads but since both show an SNP I doubt it is an artifact. 
# While more reads would be good to confirm, I believe it is a genuine SNP.


### Question 2.6 ###
#The SNP occurs between SCC2 and SAS4 and is located at chrIV:825,834.