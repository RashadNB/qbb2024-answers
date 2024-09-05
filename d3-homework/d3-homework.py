#!/usr/bin/env python3

import sys
import numpy as np

fs = open(sys.argv[1], mode='r') # Open file
fs.readline() # skip 2 lines
fs.readline()
line = (fs.readline()) 
fields = line.strip('\n').split('\t') # split column header by tabs and skip first two entries
tissues = fields[2:]
gene_names = []# create way to hold gene names
gene_IDs = [] # create way to hold gene IDs
expression = [] # create way to hold expression values

for line in fs: #for each line
    fields = line.strip('\n').split('\t') #split line
    gene_IDs.append(fields[0]) #save field 0 into gene IDs
    gene_names.append(fields[1]) #save field 1 into gene names
    expression.append(fields[2:]) #save 2+ into expression values
print()

fs.close()

# Part 2

tissues = np.array(tissues)
gene_IDs = np.array(gene_IDs)
gene_names = np.array(gene_names)
expression = np.array(expression, dtype=float)
# print(gene_IDs)
# print(gene_names)
# print(expression)

# Part 4
ten_means = np.mean(expression[0:10, :] , axis = 1)
# print(ten_means)

#Part 5

all_mean = np.mean(expression[: , :])
all_median = np.median(expression[: , :])
#print(all_mean)
#print(all_median)

#Part 6
ps_count = expression + 1
trans = np.log2(ps_count)
trans_mean = np.mean(trans[:,:])
trans_median = np.median(trans[:,:])
#print(trans_mean)
#print(trans_median)
# The median remains similar but the mean is adjusted way down to be closer to the median.

#Part 7

trans_sort = np.sort(trans, axis = 1)
c2 = trans_sort[:,-2]
c1 = trans_sort[:,-1]
c_diff = c1-c2
# print(c_diff)

#Part 8
diff_thresh = np.where(c_diff > 10)
diff_len = np.size(diff_thresh)
# print(diff_len)



