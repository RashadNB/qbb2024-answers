#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )
for my_line in my_file:
    if "##" in my_line:
        continue
    fields = my_line.split("\t")
    if fields[2] == "gene":
        fields2 = fields[8].split(";")
        gene_name = fields2[2].lstrip("gene_name \"").rstrip("\"")
        print(fields[0] + "\t" + fields[3] + "\t" + fields[4] + "\t" + gene_name)


    my_file.close