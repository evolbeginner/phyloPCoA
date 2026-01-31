iqtree2 -s cytb.aln -pre cytb -T 8 -redo -m GTR+G
iqtree2 -s cytb.aln -pre cytb -T 8 -redo -m MFP
transeq -sequence cytb.fa -outseq cytb.faa -frame 1  -table 2
mafft cytb.faa > prot.aln
iqtree2 -s prot.aln -pre prot -T 8 -redo -m MFP
iqtree2 -s prot.aln -pre prot -T 8 -redo -m LG+G 
 iqtree2 -s prot.aln -pre iqtree -T 1 -redo -m MFP -te finches_unrooted.tre
iqtree2 -s cytb.aln -pre cytb_unroot -T 8 -redo -m MFP -te finches_unrooted.tre 
