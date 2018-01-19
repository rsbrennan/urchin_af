# this script counts the number of variants falling in coding vs non-coding regions. based on
# ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3

#save variants in exons:
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 -wa -wb \
      | grep 'exon'| cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_exon.bed

# count variants in exons:
echo "number of variants in exons"
cat ~/urchin_af/analysis/variants_exon.bed | wc -l

# count variants in non-coding:
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 -wa -wb | grep -v 'exon' |\
    cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_noncoding.1.bed

# however, many of these "gene" or "transcript" annotations are exons also. So want to parse those out.

# compare file1 and file2 and then copy the rows that are present in file1 but not in file2
# this grep compares  grep -xvFf file2 file1
echo "variants in gene or transcript, but not in exon"
grep -xvFf ~/urchin_af/analysis/variants_exon.bed ~/urchin_af/analysis/variants_noncoding.1.bed | wc -l

# comm -23 <(sort file1) <(sort file2)
comm -23 <(sort ~/urchin_af/analysis/variants_noncoding.1.bed) <(sort ~/urchin_af/analysis/variants_exon.bed) | wc -l
# these two give the same number: 5867

# save output
comm -23 <(sort ~/urchin_af/analysis/variants_noncoding.1.bed) <(sort ~/urchin_af/analysis/variants_exon.bed) > ~/urchin_af/analysis/variants_noncoding.2.bed

# now find variants with no overlap in gff
zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^#' | awk '{OFS="\t"; {print "chr"$1,$2-1,$2}}' |sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 -v | > ~/urchin_af/analysis/variants_noncoding.3.bed

echo "variants in no model"
cat ~/urchin_af/analysis/variants_noncoding.3.bed | wc -l

cat ~/urchin_af/analysis/variants_noncoding.2.bed ~/urchin_af/analysis/variants_noncoding.3.bed > ~/urchin_af/analysis/variants_noncoding.bed

echo "total non-coding"

cat ~/urchin_af/analysis/variants_noncoding.bed| wc -l


#################
# selected loci
#################

# save variants in exons:
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g' | ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 -wa -wb \
      | grep 'exon'| cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_selected_exon.bed
# count variants in exons:
echo "number of significant variants in coding"
cat ~/urchin_af/analysis/variants_selected_exon.bed | wc -l

# count variants in non-coding:
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g'| ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 -wa -wb | grep -v 'exon' |\
    cut -f 1-3 | sort | uniq > ~/urchin_af/analysis/variants_selected_noncoding.1.bed

# however, many of these "gene" or "transcript" annotations are exons also. So want to parse those out.
# compare file1 and file2 and then copy the rows that are present in file1 but not in file2
# this grep compares  grep -xvFf file2 file1
echo "variants in gene or transcript, but not in exon"
grep -xvFf ~/urchin_af/analysis/variants_selected_exon.bed ~/urchin_af/analysis/variants_selected_noncoding.1.bed | wc -l

# comm -23 <(sort file1) <(sort file2)
comm -23 <(sort ~/urchin_af/analysis/variants_selected_noncoding.1.bed) <(sort ~/urchin_af/analysis/variants_selected_exon.bed) | wc -l
# these two give the same number: 112

# save output
comm -23 <(sort ~/urchin_af/analysis/variants_selected_noncoding.1.bed) <(sort ~/urchin_af/analysis/variants_selected_exon.bed) > ~/urchin_af/analysis/variants_selected_noncoding.2.bed

# now find variants with no overlap in gff
cat ~/urchin_af/analysis/cmh.selected.bed| sort -k1,1 -k2,2n | sed 's/^chr//g' | ~/bin/bedtools2/bin/bedtools intersect -a stdin -b ~/reference/SpBase3.1_build8.gff3/Transcriptome.sorted.gff3 -v | > ~/urchin_af/analysis/variants_selected_noncoding.3.bed

echo "variants in no model"
cat ~/urchin_af/analysis/variants_selected_noncoding.3.bed | wc -l

cat ~/urchin_af/analysis/variants_selected_noncoding.2.bed ~/urchin_af/analysis/variants_selected_noncoding.3.bed > ~/urchin_af/analysis/variants_selected_noncoding.bed

echo "total selected non-coding"

cat ~/urchin_af/analysis/variants_selected_noncoding.bed | wc -l




