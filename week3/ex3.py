#!/usr/bin/env python3

f = open('/Users/cmdb/qbb2024-answers/week3/biallelic.vcf') # Open the VCF file and prepare to read
count = 0 # Count the header lines that start with "##"
for line in f:
    if line.startswith("##"):
        count += 1
f.close()   

f = open('/Users/cmdb/qbb2024-answers/week3/biallelic.vcf') # Open the file again and skip the header
output = open('AF.txt', 'w')
dr = open('DP.txt', 'w')

output.write("Allele Frequency\n") # Write the headers for the output files
dr.write("Depth Reads\n")

for i, line in enumerate(f): # Skip the header lines and process the VCF content
    if i < count:
        continue  # Skip header lines
    columns = line.strip().split('\t')
    info_field = columns[7] # Access the INFO field (8th column, index 7)
    try: # Extract allele frequency from the INFO field (4th element is the AF entry)
        af_entry = info_field.split(';')[3].split('=')[1]
        output.write(af_entry + "\n")
    except IndexError:
        continue  # In case there's a problem with extracting allele frequency (I was getting an error before)
    for genotype_field in columns[9:]: # Extract depth reads from the genotype fields (9th column onwards)
        try:
            depth_read = genotype_field.split(":")[2]  # Depth is the 3rd field
            dr.write(depth_read + "\n")
        except IndexError:
            continue  # In case there's a problem with extracting depth read (I was getting an error before)

# Close all files
f.close()
output.close()
dr.close()

## Question 3.1 ##
# The output of the AF.txt file appears to be a normal (or maybe poisson) distrubtion with an average allele frequency of around 0.5.
# The appearance of the curve is expected to me since I think it would be quite uncommon for variant alleles to become fixed in any population.
# There are, however, a handful that have become fixed (frequency of 1).

## Question 3.2 ##
# Based on what we observed earlier, this distribution is expected.
# Despite the fact that we saw many unread sequences, most sit at around 4 and there are a few that get up to as high as 19 reads.
# Altogether this validates the average coverage of 4x that we obtained earlier.
# This one is definitely a poisson distribution.