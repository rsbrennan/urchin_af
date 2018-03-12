#!/bin/bash

cd ~/urchin_af/variants/

zcat ~/urchin_af/variants/urchin_filt1.vcf.gz | \
    ~/bin/vcftools/bin/vcftools --vcf - \
    --positions ~/urchin_af/variants/keep.snps.txt \
    --remove-indv OASV2_DNA_D1_7_5_S_07 \
    --remove-indv OASV2_DNA_D1_7_5_S_03 \
    --remove-indv OASV2_DNA_D1_7_5_S_04 \
    --remove-indv OASV2_DNA_D1_7_5_S_05 \
    --recode -c | \
    bgzip > ~/urchin_af/variants/urchin_filt2.vcf.gz

tabix -p vcf ~/urchin_af/variants/urchin_filt2.vcf.gz
