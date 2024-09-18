#!/usr/bin/env python3


import sys # import packages
import numpy


reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT'] # generate list of sequence reads.


k = 3 # each k-mer has 3 nucleotides so setting k = 3.


graph = [] # create empty list.

# generate edge nodes of de Bruijn graph and add them to the empty "graph" list
for read in reads:
    for i in range(len(read) - k):                          # iterate through the list at indices 0 and 1...
        kmer1 = read[i : i + k]                             # kmer1 is substring of nucleotides starting from position i...
        kmer2 = read[i + 1 : i + 1 + k]                     # kmer2 is the substring of nucleotides starting from position i+1...
        graph.append(kmer1 + " -> " + kmer2)                # add overlapping 3-mers of kmer1 and kmer2 to the list as edge nodes.
        
# print(graph)



header = "digraph {"      # open graph description.
print(header)             

for edge in graph:        # print each edge of the graph.
    print(edge)                                             


footer = "}"              # close graph description
print(footer)
