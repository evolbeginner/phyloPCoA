le metadata1.txt |perl -e ' $i =0; while(<>) {chomp; my @a = split/\s+/;$i=$i+1; print "$i\t$a[0]\t$a[1]_$a[2]\n"; } ' |le > metadata.txt
