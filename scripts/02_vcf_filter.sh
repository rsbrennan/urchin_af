#!/bin/bash

cd ~/urchin_af/variants/

# concatenate all output
    # use vcftools, vcf-concat

    #do the following bc error in location of perl library. this fixes it.
    export PERL5LIB=$PERL5LIB:~/bin/vcftools/src/perl

    ~/bin/vcftools/bin/vcf-concat $(ls -1 ~/urchin_af/variants/temp/*.vcf | perl -pe 's/\n/ /g')> ~/urchin_af/variants/tmp.vcf

        ~/bin/vcftools/bin/vcftools --vcf ~/urchin_af/variants/tmp.vcf \
        --maf 0.01 \
        --minQ 20 \
        --min-alleles 2 --max-alleles 2 \
        --remove-indels --min-meanDP 10 \
        --max-missing .9 \
        --recode -c | \
        ~/bin/vcftools/bin/vcf-sort -c | \
        bgzip >  ~/urchin_af/variants/urchin_dupsincl.vcf.gz

    tabix -p vcf ~/urchin_af/variants/urchin_dupsincl.vcf.gz



