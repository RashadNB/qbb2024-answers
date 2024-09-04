#!/usr/bin/env python3


import sys

my_file = open(sys.argv[2])
for line in my_file:
    line = line.rstrip('\n')
    if sys.argv[1] in line:
        print(line)
# ./grep.py FIS1 gencode.v46.basic.annotation.gtf |less -S 
# ^ unix code corresponding to python code