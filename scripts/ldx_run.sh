
cd ~/urchin_af/variants/

samtools=~/bin/samtools-1.6/samtools
vcf_i=~/urchin_af/variants/urchin_final.vcf
#######################
# loop over samples
#######################
for rep in `$samtools view -H /data/OASV2/bwamem/fullymapped/bwamem_final.merged.sorted.bam | grep "^@RG" | cut -f 2 | cut -f 2 -d ":" | cut -f 1 -d "-" | sort | uniq`; do

echo $rep

bam_in=~/urchin_af/variants/${rep}.bam

# check if outfile exists, if it does, delete.
if [ -e ~/urchin_af/analysis/${rep}.ld.out ]
then
    echo "${rep}.ld.out exists, removing now"
    rm ~/urchin_af/analysis/${rep}.ld.out
else
    echo "${rep}.ld.out does not exist, doing nothing"

fi

#cycling over each scaffold, to be memory efficient.

cat $vcf_i | grep -v '^#'| cut -f 1 | sort | uniq > ~/urchin_af/variants/scaffolds.txt

cat $vcf_i | grep '^#' > header.txt

while read scaf; do

        #subset to correct scaffold
        $samtools view -h $bam_in $scaf > tmp.${rep}.scaf.sam

        #subset vcf
        cat $vcf_i | grep -v '^#' | grep -w $scaf > tmp1.${rep}.vcf
        cat  header.txt tmp1.${rep}.vcf > tmp.${rep}.vcf

        #calculate ld on the scaffold
        ~/bin/LDx.pl -l 40 -h 4800 -s 500 -q 20 -a 0.05 -i 5 tmp.${rep}.scaf.sam tmp.${rep}.vcf > ~/urchin_af/analysis/ld.${rep}.tmp

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
