#!/usr/bin/env python3

#Import packages
import sys
import numpy as np

filename = sys.argv[1] #get gene-tisse file name. Opens file that is the first argument of command line.
fs = open(filename, mode='r')  #open gene-tissue file.
relevant_samples = {} #create dictionary to hold samples for gene-tissue pairs.
for line in fs: #step through file.
    fields = line.rstrip('\n').split('\t') #Split line in file on tabs.
    key = (fields[0], fields[2]) #create key from gene and tissue.
    relevant_samples[key] = [] #initialize dictionary from key with list to hold samples.
fs.close() #close the file.

#print(relevant_samples)

filename = sys.argv[2] #get metadata file name
fs = open(filename, mode='r')  #open file
fs.readline() #skipping the first line
tissue_samples = {} #create dictionary to hold samples for gene-tissue pairs
for line in fs: #step through file
    fields = line.rstrip('\n').split('\t') #Split line into fields
    key = fields[6] #assign keys based on 6th field (SMTSD)
    value = fields[0] #assign values based on the 0th field (Sample IDs)
    tissue_samples.setdefault(key, []) #if you hit a tissue you haven't seen before, add it to the dictionary with the value of a blank list.
    tissue_samples[key].append(value) #add tissues you've seen before to the list assigned to them in the previous line.
fs.close() #close the file

#print(tissue_samples) #print the dictionary "tissue_samples"

filename = sys.argv[3] #get RNASeQ file name
fs = open(filename, mode='r')  #open file
fs.readline() #skipping the first line
fs.readline() #skipping the second line
header = fs.readline().rstrip('\n').split('\t') #just taking the 3rd line because we read the first two in the previous two lines. Only does it once because it's not part of a for loop.
header = header[2:] # getting sample IDs.
#print(header)

tissue_columns = {} #making a dictionary to hold indices of samples associated with tissues.
for tissue, samples in tissue_samples.items(): #items retreives both values and keys from the dictionary we previously made.
    tissue_columns.setdefault(tissue, []) #every time we encounter a new tissue we make a new list for it.
    for sample in samples: #for every sample value that we retrieved above...
        if sample in header: #if it's in the document header...
            position = header.index(sample) #record it's index value in the document...
            tissue_columns[tissue].append(position) #and add it to the dictionary.
#print(tissue_columns) #print the dictionary "tissue columns".

#for key, value in tissue_columns.items():#get the length of the list of sample columns.
    #print(key, len(value)) #printing the key to each sample list along with the length of each list.

max_value = 0 #Find tissue with max number of samples
max_tissue = ""

for tissue, samples in tissue_columns.items():
    if len(samples) >= max_value:
        max_value = len(samples)
        max_tissue = tissue
#print(max_tissue, max_value)

min_value = 0
min_tissue = ""
for tissue, samples in tissue_columns.items():
    if len(samples) <= min_value:
        min_value = len(samples)
        min_tissue = tissue
print(min_tissue, min_value)
# Whole blood, frontal brain cortex, and subcutaneous adipose tissue have the largest number of samples.
# Bladder, kidney medulla, and CML have the least samples.





