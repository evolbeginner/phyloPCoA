le metadata1.txt |perl -e ' $i =0; while(<>) {chomp; my @a = split/\t/; if ($a[2] eq "H" ){ $i=$i+1; print "$i\t$a[0]\t$a[1]\n";} } ' > metadata_H.txt
le metadata1.txt |perl -e ' $i =0; while(<>) {chomp; my @a = split/\t/; if ($a[2] eq "F" ){ $i=$i+1; print "$i\t$a[0]\t$a[1]\n";} } ' > metadata_F.txt
le metadata1.txt |perl -e ' $i =0; while(<>) {chomp; my @a = split/\t/; if ($a[2] eq "M" ){ $i=$i+1; print "$i\t$a[0]\t$a[1]\n";} } ' > metadata_M.txt
python split_abundance_by_group.py 
