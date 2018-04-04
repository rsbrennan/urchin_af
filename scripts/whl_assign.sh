# make bedfile from cmh file: ~/urchin_af/analysis/cmh.master.out

cat ~/urchin_af/analysis/cmh.master.out | tail -n +2 | awk '{OFS="\t"} { print $1, $2-1, $2 }' | sort -k1,1 -k2,2n > ~/urchin_af/analysis/cmh.all.sorted.bed

cat ~/urchin_af/analysis/cmh.master.out  |sort -k1,1 -k2,2n > ~/urchin_af/analysis/cmh.master.sort.out


### assign whl names to genes

# all cmh
~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.all.sorted.bed \
-b ~/reference/urchin_probes.sorted.bed \
-d -t first  > ~/urchin_af/analysis/cmh.all.genes.txt

# convert cmh.genes.txt gene names to those compatible with echinobase
#get id from 5th column

# remove transcript number from whl names
cat ~/urchin_af/analysis/cmh.all.genes.txt | grep -ohP "ID=.*?;" | sed 's/ID=//g' | sed 's/;//g' | sed 's/"//g'| \
    cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.all.id

