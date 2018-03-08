#!/bin/bash/ -l

cd ~/urchin_af/variants/


#generate bedfile for each probe + 600bp in either direction:
cat /data/OASV2/RG_7601_Probe.extended.bed | awk '{OFS="\t"; if (!/^#/){print $1,$2-600,$3+600}}' | awk -v OFS='\t' '$2<0 {$2=0} 1' > ~/urchin_af/analysis/probe_pi_keep.bed

#then bed file for selected +- 600

 sort -k1,1 -k2,2n ~/urchin_af/analysis/cmh.selected.bed |   bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.bed > ~/urchin_af/analysis/selected_probe_keep.bed

for p in D1_7.pop D7_7.pop D1_8.pop D7_8.pop; do

    ~/bin/vcftools/bin/vcftools --gzvcf ~/urchin_af/variants/urchin_final.vcf.gz  \
    --keep $p --bed ~/urchin_af/analysis/probe_pi_keep.bed \
    --window-pi 400 --window-pi-step 200  \
    --out ~/urchin_af/analysis/${p}

    ~/bin/vcftools/bin/vcftools --gzvcf ~/urchin_af/variants/urchin_final.vcf.gz  \
    --keep $p --bed ~/urchin_af/analysis/selected_probe_keep.bed \
    --window-pi 400 --window-pi-step 200 \
    --out ~/urchin_af/analysis/${p}.selected

    echo $p

done
