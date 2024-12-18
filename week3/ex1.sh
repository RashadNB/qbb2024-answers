#!/usr/bin/env bash

#!/usr/bin/env bash

### Question 1.1 ###
wget https://www.dropbox.com/s/ogx7nbhlhrp3ul6/BYxRM.tar.gz # download dropbox file.
tar -xvzf BYxRM.tar.gz # unzip file.
FastQC A01_09.fastq
# the sequencing reads are 76 bp long.


### Question 1.2 ###
# Divide by four because each read takes up 4 lines.
wc -l A01_09.fastq # 2678192
expr 2678192 / 4  # 669548
# There are 669548 reads in the file.


### Question 1.3 ###
# find size of reference genome by grepping for lines that start with > which are the sequences.
grep -v '^>' sacCer3.fa | tr -d '\n' | wc -c

# reference genome is 12157105 bp long.

# expected average depth of coverage genome length divided by reads
echo "76 * 669548 / 12200000" | bc 
# the depth of coverage is 4

### Question 1.4 ###
du -m A01_* | sort -n
# A01_62.fastq file is the biggest (149Mb du), A01_27.fastq is the smallest (110Mb of du).

### Question 1.5 ###
FastQC *.fastq

# The median base quality on the read is 37.
# This means that the probability of any base being an error is 10e-3.7.
# The variability of quality seems to increase at both ends of a given read.