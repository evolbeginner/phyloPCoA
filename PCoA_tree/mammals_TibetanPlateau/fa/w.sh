transeq -sequence cytb.fa -outseq cytb.faa -frame 1  -table 2
mafft cytb.faa > cytb.aln
iqtree2 -s cytb.aln  -pre cytb -T 8 -redo -m LG+G
iqtree2 -s cytb.aln -pre cytb -T 8 -redo -m MFP
 iqtree2 -s cytb.aln -pre iqtree -T 1 -redo -m MFP -te mammals_unrooted.tre 
