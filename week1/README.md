#Answer 1
- 1 Mbp * 3 = 3 million / 100 = 30,000 reads are required.

#Answer 1.4
- Approximately 50,000 sites have not been sequenced.
- Both the normal and poisson distrubutions are not perfect, but they model the data fairly well. They fall a little short of the top of the curve. The curves corroborate each other.

#Answer 1.5
- grep-ing to remove 10s, 20s, and 30s, and then grepping for 0s shows that 60 sites were not covered.
- 'grep -v 10 coverage2.txt | grep -v 20 | grep -v 30 |  grep 0 | wc -l'
- Both the normal and poisson distrubutions are not perfect, but they model the data fairly well. They fall a little short of the top of the curve. The curves corroborate each other.

#Answer 1.6
- Doing the same as above but including 40s, 50s, and 60s shows that there are 3 uncovered sites.
- 'grep -v 10 coverage3.txt | grep -v 20 | grep -v 30 | grep -v 40 | grep -v 50 | grep -v 60 | grep -v 70 | grep 0 | wc -l'
- Once again both distributions match the data quite well.

#Answer 2.4
- 'dot -Tpng deBruijn_graph.dot -o deBruijn_graph.png' 

#Answer 2.5
- 5' TCTTATTCATTGATTT 3'

#Answer 2.6
- The main thing needed to accurately reconstruct the genome would be longer sequence reads. Longer reads would allow for more overlap between fragments rather than the small number of nucleotide overlaps we can get with our current set-up. This would reduce the number of possible subsequent sequences. 