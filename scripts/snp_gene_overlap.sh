
#sort and pull out id
sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.all_dups_in.genes.txt > ~/urchin_af/analysis/cmh.all_dups_in.sorted.genes.txt

cat ~/urchin_af/analysis/cmh.all_dups_in.sorted.genes.txt | cut -f 14 > ~/urchin_af/analysis/cmh.all_dups_in.ID.txt

sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.selected.bed > ~/urchin_af/analysis/cmh.selected.sorted.bed
sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.neutral.bed > ~/urchin_af/analysis/cmh.neutral.sorted.bed
sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.all.bed > ~/urchin_af/analysis/cmh.all.sorted.bed

~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.selected.sorted.bed \
-b ~/reference/urchin_probes.sorted.bed \
-d -t first  > ~/urchin_af/analysis/cmh.selected.genes.txt

~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.neutral.sorted.bed \
-b ~/reference/urchin_probes.sorted.bed \
-d -t first  > ~/urchin_af/analysis/cmh.neutral.genes.txt

# all cmh
~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.all.sorted.bed \
-b ~/reference/urchin_probes.sorted.bed \
-d -t first  > ~/urchin_af/analysis/cmh.all.genes.txt

# convert cmh.genes.txt gene names to those compatible with echinobase
#get id from 5th column

# remove transcript number from whl names
#some snps fall in more than one transcript- it looks like different versions of the same gene. In this case, just using the base transcript.

cat ~/urchin_af/analysis/cmh.genes.txt | grep -ohP "ID=.*?;" | sed 's/ID=//g' | sed 's/;//g' | sed 's/"//g'| \
    cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.selected.id

cat ~/urchin_af/analysis/cmh.neutral.genes.txt | grep -ohP "ID=.*?;" | sed 's/ID=//g' | sed 's/;//g' | sed 's/"//g'| \
    cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.neutral.id

cat ~/urchin_af/analysis/cmh.all.genes.txt | grep -ohP "ID=.*?;" | sed 's/ID=//g' | sed 's/;//g' | sed 's/"//g'| \
    cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.all.id
