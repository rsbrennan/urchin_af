#!/bin/bash

cd ~/urchin_af/variants/

zcat ~/urchin_af/variants/urchin_dupsincl.vcf.gz | \
    ~/bin/vcftools/bin/vcftools --vcf - \
    --positions ~/urchin_af/variants/keep.snps.txt \
    --remove-indv OASV2_DNA_D1_7_5_S_07 \
    --recode -c | \
    bgzip > ~/urchin_af/variants/urchin_filtered.vcf.gz

tabix -p vcf ~/urchin_af/variants/urchin_filtered.vcf.gz
