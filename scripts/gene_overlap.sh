#!/bin/bash

sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.selected.bed > ~/urchin_af/analysis/cmh.selected.sorted.bed
sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.neutral.bed > ~/urchin_af/analysis/cmh.neutral.sorted.bed
sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.all.bed > ~/urchin_af/analysis/cmh.all.sorted.bed
sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.out.txt > ~/urchin_af/analysis/cmh.out.sorted.txt
sort -k1,1 -k4,4n ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 > ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.1.gff3

#selected
~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.selected.sorted.bed \
-b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.1.gff3 \
-d -t all > ~/urchin_af/analysis/cmh.selected.genes.txt

#neutral genes
~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.neutral.sorted.bed \
-b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.1.gff3 \
-d -t all > ~/urchin_af/analysis/cmh.neutral.genes.txt

# all cmh
~/bin/bedtools2/bin/bedtools closest -a ~/urchin_af/analysis/cmh.all.sorted.bed \
-b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.1.gff3 \
-d -t all > ~/urchin_af/analysis/cmh.all.genes.txt

###
###### make ID column from gene files. Use this file for GO and SPU assignment
###

# isolate whl ids
cat ~/urchin_af/analysis/cmh.selected.genes.txt | cut -f 12 | cut -f 1 -d ";" | sed 's/ID=//g' | sed 's/"//g' | cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.selected.id

cat ~/urchin_af/analysis/cmh.neutral.genes.txt | cut -f 12 | cut -f 1 -d ";" | sed 's/ID=//g' | sed 's/"//g' | cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.neutral.id

cat ~/urchin_af/analysis/cmh.all.genes.txt | cut -f 12 | cut -f 1 -d ";" | sed 's/ID=//g' | sed 's/"//g' | cut -f 1-2 -d '.' > ~/urchin_af/analysis/cmh.all.id

paste  ~/urchin_af/analysis/cmh.selected.genes.txt ~/urchin_af/analysis/cmh.selected.id > ~/urchin_af/analysis/cmh.selected.genes.2.txt
paste  ~/urchin_af/analysis/cmh.neutral.genes.txt ~/urchin_af/analysis/cmh.neutral.id > ~/urchin_af/analysis/cmh.neutral.genes.2.txt
paste  ~/urchin_af/analysis/cmh.all.genes.txt ~/urchin_af/analysis/cmh.all.id > ~/urchin_af/analysis/cmh.all.genes.2.txt

mv ~/urchin_af/analysis/cmh.selected.genes.2.txt ~/urchin_af/analysis/cmh.selected.genes.txt
mv ~/urchin_af/analysis/cmh.neutral.genes.2.txt ~/urchin_af/analysis/cmh.neutral.genes.txt
mv ~/urchin_af/analysis/cmh.all.genes.2.txt ~/urchin_af/analysis/cmh.all.genes.txt
