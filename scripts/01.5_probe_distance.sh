#!/bin/bash

# sort the probe bed file
# this file contains all probes used to capture dna and make libraries
sort -k1,1 -k2,2n /data/OASV2/RG_7601_Probe.extended.bed > ~/reference/urchin_probes.sorted.bed

#bedfile from vcf
zcat ~/urchin_af/variants/urchin_filt2.vcf.gz | grep -v '^#' | cut -f 1-2 |\
awk '{OFS="\t";{print $1,$2-1,$2}}' |\
tail -n +2 | sort -k1,1 -k2,2n > ~/urchin_af/data/variant.probedist.bed

# calc distance
bedtools closest -a ~/urchin_af/data/variant.probedist.bed \
-b ~/reference/urchin_probes.sorted.bed \
-d -t first  > ~/urchin_af/analysis/snp_probe_dist.txt
