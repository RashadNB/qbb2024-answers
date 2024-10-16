#!/usr/bin/env python3

#grep.py     target      myFile.txt
#sys.argv[0]  sys.argv[1]  sys.argv[2]

#grep finds strings in a file. It will give you the whole line in which your string is found unless you remo

import sys #sys is a package that allows you to grep

my_file = open(sys.argv[2]) #this is a list of all the arguments in the command line. In this case it is OPENING the 2nd item in the terminal which is the file genes.gtf
for line in my_file:
    line = line.rstrip('\n')
    if sys.argv[1] in line:
        print(line)
# ./grep.py FIS1 gencode.v46.basic.annotation.gtf |less -S 
# ^ unix code corresponding to python code