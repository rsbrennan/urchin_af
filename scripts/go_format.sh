

for snp_set in master genic exon intron intergenic non_coding; do

    #these are all temporary and can be written over and removed
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out | cut -f 3 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/snp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out | cut -f 3 \
         > ~/urchin_af/analysis/go_enrichment/snp_head.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out | cut -f 14 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/GO.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out | cut -f 6 \
        > ~/urchin_af/analysis/go_enrichment/sig.out

    cat ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out | cut -f 15  \
        > ~/urchin_af/analysis/go_enrichment/logp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out | cut -f 4  \
        > ~/urchin_af/analysis/go_enrichment/pval.out

    # make go annotation with no header
    paste ~/urchin_af/analysis/go_enrichment/snp.out ~/urchin_af/analysis/go_enrichment/cmh.pH75_${snp_set}_GO.out \
        > ~/urchin_af/analysis/go_enrichment/GO.pH75_${snp_set}.annotation

    # make significance tables with headers
    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/sig.out -d "," | \
         tr -d " \t"    | sed 's/FALSE/0/g' | sed 's/TRUE/1/g' \
         > ~/urchin_af/analysis/go_enrichment/sig_binary.pH75.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/logp.out -d "," | \
        tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/logp.pH75.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/pval.out -d "," | \
         tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/pval.pH75.${snp_set}.table

    # change for topGO format, need to be comma delimted
    cat ~/urchin_af/analysis/go_enrichment/GO.${snp_set}.annotation | sed 's/;/,/g' | \
        grep -v 'unknown' > ~/urchin_af/analysis/go_enrichment/topGO.pH75.${snp_set}.annotation

done


for snp_set in master genic exon intron intergenic non_coding; do

    #these are all temporary and can be written over and removed
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out | cut -f 3 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/snp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out | cut -f 3 \
         > ~/urchin_af/analysis/go_enrichment/snp_head.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out | cut -f 14 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/GO.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out | cut -f 6 \
        > ~/urchin_af/analysis/go_enrichment/sig.out

    cat ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out | cut -f 15  \
        > ~/urchin_af/analysis/go_enrichment/logp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out | cut -f 4  \
        > ~/urchin_af/analysis/go_enrichment/pval.out

    # make go annotation with no header
    paste ~/urchin_af/analysis/go_enrichment/snp.out ~/urchin_af/analysis/go_enrichment/cmh.pH80_${snp_set}_GO.out \
        > ~/urchin_af/analysis/go_enrichment/GO.pH80_${snp_set}.annotation

    # make significance tables with headers
    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/sig.out -d "," | \
         tr -d " \t"    | sed 's/FALSE/0/g' | sed 's/TRUE/1/g' \
         > ~/urchin_af/analysis/go_enrichment/sig_binary.pH80.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/logp.out -d "," | \
        tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/logp.pH80.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/pval.out -d "," | \
         tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/pval.pH80.${snp_set}.table

    # change for topGO format, need to be comma delimted
    cat ~/urchin_af/analysis/go_enrichment/GO.${snp_set}.annotation | sed 's/;/,/g' | \
        grep -v 'unknown' > ~/urchin_af/analysis/go_enrichment/topGO.pH80.${snp_set}.annotation

done

for snp_set in master; do

    #these are all temporary and can be written over and removed
    cat ~/urchin_af/analysis/go_enrichment/cmh.overlap_${snp_set}_GO.out | cut -f 3 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/snp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.overlap_${snp_set}_GO.out | cut -f 3 \
         > ~/urchin_af/analysis/go_enrichment/snp_head.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.overlap_${snp_set}_GO.out | cut -f 14 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/GO.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.overlap_${snp_set}_GO.out | cut -f 6 \
        > ~/urchin_af/analysis/go_enrichment/sig.out

    cat ~/urchin_af/analysis/go_enrichment/cmh.overlap_${snp_set}_GO.out | cut -f 15  \
        > ~/urchin_af/analysis/go_enrichment/logp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.overlap_${snp_set}_GO.out | cut -f 4  \
        > ~/urchin_af/analysis/go_enrichment/pval.out

    # make go annotation with no header
    paste ~/urchin_af/analysis/go_enrichment/snp.out ~/urchin_af/analysis/go_enrichment/GO.out \
        > ~/urchin_af/analysis/go_enrichment/GO.overlap_${snp_set}.annotation

    # make significance tables with headers
    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/sig.out -d "," | \
         tr -d " \t"    | sed 's/FALSE/0/g' | sed 's/TRUE/1/g' \
         > ~/urchin_af/analysis/go_enrichment/sig_binary.overlap.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/logp.out -d "," | \
        tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/logp.overlap.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/pval.out -d "," | \
         tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/pval.overlap.${snp_set}.table

    # change for topGO format, need to be comma delimted
    cat ~/urchin_af/analysis/go_enrichment/GO.overlap_${snp_set}.annotation | sed 's/;/,/g' | \
         > ~/urchin_af/analysis/go_enrichment/topGO.overlap.${snp_set}.annotation

done
