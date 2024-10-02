#!/usr/bin/env bash 

### Question 2.1 ###
grep -c "^>" sacCer3.fa # grep for > because chromosomes start with > in FASTA files.
# There are 17 chromosomes in the yeast genome (16 plus one mitochondrial).

