#!/bin/bash/ -l

cd ~/urchin_af/variants/

for p in D1_7.pop D7_7.pop D1_8.pop D7_8.pop; do

    ~/bin/vcftools/bin/vcftools --gzvcf ~/urchin_af/variants/urchin_final.vcf.gz  \
    --keep $p --bed ~/urchin_af/analysis/probe_pi_keep.bed \
    --window-pi 150 --window-pi-step 50  \
    --out ~/urchin_af/analysis/${p}

    ~/bin/vcftools/bin/vcftools --gzvcf ~/urchin_af/variants/urchin_final.vcf.gz  \
    --keep $p --bed ~/urchin_af/analysis/selected_probe_keep.bed \
    --window-pi 150 --window-pi-step 50 \
    --out ~/urchin_af/analysis/${p}.selected

    echo $p

done
