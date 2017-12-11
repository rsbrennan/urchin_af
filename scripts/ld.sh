#!/bin/bash/ -l

#using LDx to calculate LD

samtools=~/bin/samtools-1.6/samtools

# generate sam files for analysis

cd ~/urchin_af/variants/

$samtools view -H /data/OASV2/merged.fixmate.sorted.bam > ~/urchin_af/variants/sam.header

#check for presence of D1 bam
if [[ -s ~/urchin_af/variants/D1.bam ]];
then
    echo "D1.bam exists and is not empty. do nothing"
else
    echo "D1.bam doesn't exist; creating now"
    $samtools view /data/OASV2/merged.fixmate.sorted.bam   | grep 'OASV2_DNA_D1_8_0.*' | cat ~/urchin_af/variants/sam.header - | $samtools view -Sb > ~/urchin_af/variants/D1.bam
    $samtools index ~/urchin_af/variants/D1.bam
    echo "D1 bam created"
fi

#check for presence of D7_8 bam
if [ -s ~/urchin_af/variants/D7_8.bam ]
then
     echo "D7_8.bam exists and is not empty. do nothing"
else
    echo "D7_8.bam doesn't exist; creating now"
    $samtools view /data/OASV2/merged.fixmate.sorted.bam   | grep 'OASV2_DNA_D7_8_0.*' | cat ~/urchin_af/variants/sam.header - | $samtools view -Sb > ~/urchin_af/variants/D7_8.bam
    $samtools index ~/urchin_af/variants/D7_8.bam
    echo "D7_8 bam created"
fi


#check for presence of D8_7 bam
if [ -s ~/urchin_af/variants/D7_7.bam ]
then
    echo "D7_7.bam exists and is not empty. do nothing"
else
    echo "D7_7.bam doesn't exist; creating now"
    $samtools view /data/OASV2/merged.fixmate.sorted.bam   | grep 'OASV2_DNA_D7_7_5.*' | cat ~/urchin_af/variants/sam.header - | $samtools view -Sb > ~/urchin_af/variants/D7_7.bam
    $samtools index ~/urchin_af/variants/D7_7.bam
    echo "D7_7 bam created"
fi


cd ~/urchin_af/analysis/

# input.sam - Sorted, duplicates removed .sam file of a single chromosome
# snplist.vcf - VCF file of probable SNP locations. snplist.vcf MUST be sorted.

for rep in D1 D7_7 D7_8; do

echo $rep

# check if vcf exists, if not, create.
if [ -e ~/urchin_af/variants/urchin_final.vcf ]
then
    echo "vcf exists"
else
    echo "vcf does not exist, making now"

    zcat ~/urchin_af/variants/urchin_final.vcf.gz > ~/urchin_af/variants/urchin_filtered.vcf

fi

vcf_i=~/urchin_af/variants/urchin_filtered.vcf
bam_in=~/urchin_af/variants/${rep}.bam

# check if outfile exists, if it does, delete.
if [ -e ~/urchin_af/analysis/${rep}.ld.out ]
then
    echo "${rep}.ld.out exists, removing now"
    rm ~/urchin_af/analysis/${rep}.ld.out
else
    echo "${rep}.ld.out does not exist, doing nothing"

fi


cat $vcf_i | grep -v '^#'| cut -f 1 | sort | uniq > ~/urchin_af/variants/scaffolds.txt

cat $vcf_i | grep '^#' > header.txt

while read scaf; do

        #echo $scaf

        #subset to correct scaffold
        $samtools view -h $bam_in $scaf > tmp.${rep}.scaf.sam

        #subset vcf
        cat $vcf_i | grep -v '^#' | grep -w $scaf > tmp1.${rep}.vcf
        cat  header.txt tmp1.${rep}.vcf > tmp.${rep}.vcf

        #calculate ld on the scaffold
        ~/bin/LDx.pl -l 40 -h 4800 -s 600 -q 10 -a 0.05 -i 5 tmp.${rep}.scaf.sam tmp.${rep}.vcf > ~/urchin_af/analysis/ld.${rep}.tmp

        # add scaffold identifier

        sed -i "s/^/$scaf\t/g" ~/urchin_af/analysis/ld.${rep}.tmp

        cat  ~/urchin_af/analysis/ld.${rep}.tmp >> ~/urchin_af/analysis/${rep}.ld.out

done<~/urchin_af/variants/scaffolds.txt

rm ld.${rep}.tmp
rm tmp.${rep}.vcf
rm tmp.${rep}.scaf.sam

echo $rep
echo "done"

done
