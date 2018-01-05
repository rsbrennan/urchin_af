cat ~/reference/urchin_probes.sorted.bed | awk 'BEGIN{FS=OFS="\t"} {gsub(" ", "", $1)}1' > ~/reference/urchin_probes.sorted.nospace.bed

echo "number of variants in exons"
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="exon" {print $0}'| wc -l

echo "number of exon probes represented"
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="exon" {print $0}' | cut -f 7 | sort | uniq | wc -l

echo "number of variants in regulatory"
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="regulatory" {print $0}'| wc -l

echo "number of regulatory probes represented"
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="regulatory" {print $0}'| cut -f 7 | sort | uniq | wc -l

zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first



######
##
## look at prop of sig genes in coding/non
##
#######

echo "number of significant variants in regulatory"
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g'|  ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="regulatory" {print $0}' | wc -l

echo "number of regulatory probes represented for sig variants"
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g'|  ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="regulatory" {print $0}' | cut -f 7 | sort | uniq | wc -l

echo "number of significant variants in regulatory"
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g'|  ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="exon" {print $0}' | wc -l

echo "number of regulatory probes represented for sig variants"
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g'|  ~/bin/bedtools2/bin/bedtools closest -a stdin -b ~/reference/urchin_probes.sorted.nospace.bed \
    -d -t first | awk ' $9=="exon" {print $0}' | cut -f 7 | sort | uniq | wc -l
