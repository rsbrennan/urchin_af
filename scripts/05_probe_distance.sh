#!/bin/bash

sort -k1,1 -k2,2n /data/OASV2/RG_7601_Probe.extended.bed > ~/reference/urchin_probes.sorted.bed
sort -k1,1 -k2,2n ~/reference/SpBase3.1_build8.gff3/Transcriptome.gff3 > ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3

#remove EC lines
#bedfile from AF file
awk '{OFS="\t";{print $1,$2-1,$2}}' ~/urchin_af/data/allele.freq_all.txt |\
tail -n +2 | sort -k1,1 -k2,2n > ~/urchin_af/data/allele.freq.bed


zcat ~/reference/blast2go-whl.annot.txt.gz | grep -v 'EC' > ~/reference/blast2go-whl.annot.txt.1.gz
bgzip -c ~/reference/blast2go-whl.annot.txt.1.gz > ~/reference/blast2go-whl.annot.txt.gz
rm ~/reference/blast2go-whl.annot.txt.1.gz

# calc distance
bedtools closest -a ~/urchin_af/data/allele.freq.bed \
-b ~/reference/urchin_probes.sorted.bed \
-d -t first  > ~/urchin_af/analysis/snp_probe_dist.txt
