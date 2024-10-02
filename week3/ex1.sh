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


