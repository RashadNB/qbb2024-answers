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

