# this script counts the number of variants falling in coding vs non-coding regions. based on
# the Strongylocentrotus_purpuratus.Spur_3.1.37.gff3.gz gff.

zcat ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37.gff3.gz  | grep -v 'supercontig' | sort -k1,1 -k4,4n > ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3

# count variants in exons:
echo "number of variants coding"

zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -wa -wb \
      | grep -v 'three_prime_UTR' | cut -f 1-3 | sort | uniq | wc -l

# save variants in exons:
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -wa -wb \
      | cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_exon.bed

# count variants in non-coding:
echo "number of variants non-coding"
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -v | wc -l

# save variants in non-coding:
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -v > ~/urchin_af/analysis/variants_noncoding.bed

###
# selected loci
##

# count selected in exons
echo "number of significant variants in coding"
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g' | ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -wa -wb \
    | cut -f 1-3 | sort | uniq | wc -l

cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g' | ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -wa -wb \
    | cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_selected_exon.bed

echo "number of significant variants in regulatory"
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g' | ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -v \
    | cut -f 1-3 | sort | uniq | wc -l

cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g' | ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/Strongylocentrotus_purpuratus.Spur_3.1.37_sorted.gff3 -v \
    | cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_selected_noncoding.bed
