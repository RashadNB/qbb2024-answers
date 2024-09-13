#!/usr/bin/env python3

import sys
import numpy as np
import scipy as sp

genome_size = 1000000   # 1 Mbp genome
read_length = 100  
coverage_depth = int(sys.argv[1])   #sys.argv allows me to modify my coverage depth in the terminal without altering code.

total_bases_to_sequence = genome_size * coverage_depth
number_of_reads = total_bases_to_sequence // read_length

coverage = np.zeros(genome_size, dtype=int) # Initialize the genome coverage array

start_positions = np.random.randint(0, genome_size - read_length + 1)

for _ in range(number_of_reads): # For loop to randomly select start positions, calculate end positions, and increment coverage
    start = np.random.randint(0, genome_size - read_length + 1) # Randomly select a start position (uniformly between 0 and genome_size - read_length)
    end = start + read_length # Calculate the end position (end position is start + read_length)
    coverage[start:end] +=1

for bp in coverage: # for every base pair in coverage, print the base pair
    print(bp)
