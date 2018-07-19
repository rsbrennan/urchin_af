#!/bin/bash


# filter by max depth calculated in 01.3_filter.R
cd ~/urchin_af/variants/

zcat ~/urchin_af/variants/urchin_filt1.vcf.gz | \
    ~/bin/vcftools/bin/vcftools --vcf - \
    --max-meanDP 372 \
    --recode -c | \
    bgzip > ~/urchin_af/variants/urchin_filt2.vcf.gz

tabix -p vcf ~/urchin_af/variants/urchin_filt2.vcf.gz
