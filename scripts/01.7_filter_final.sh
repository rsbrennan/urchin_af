#!/bin/bash


#removing all variants that are > 2kb from a probe.

cd ~/urchin_af/variants/

zcat ~/urchin_af/variants/urchin_filt2.vcf.gz | \
    ~/bin/vcftools/bin/vcftools --vcf - \
    --positions ~/urchin_af/variants/keep.ontarget.txt \
    --recode -c | \
    bgzip > ~/urchin_af/variants/urchin_final.vcf.gz

tabix -p vcf ~/urchin_af/variants/urchin_final.vcf.gz
