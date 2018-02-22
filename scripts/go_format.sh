for snp_set in master genic exon synonymous non_syn intron intergenic non_coding; do

    #these are all temporary and can be written over and removed
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}_GO.out | cut -f 3 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/snp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}_GO.out | cut -f 3 \
         > ~/urchin_af/analysis/go_enrichment/snp_head.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}_GO.out | cut -f 12 | \
        tail -n +2 > ~/urchin_af/analysis/go_enrichment/GO.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}_GO.out | cut -f 5 \
        > ~/urchin_af/analysis/go_enrichment/sig.out

    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}_GO.out | cut -f 13  \
        > ~/urchin_af/analysis/go_enrichment/logp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}_GO.out | cut -f 4  \
        > ~/urchin_af/analysis/go_enrichment/pval.out

    # make go annotation with no header
    paste ~/urchin_af/analysis/go_enrichment/snp.out ~/urchin_af/analysis/go_enrichment/GO.out \
        > ~/urchin_af/analysis/go_enrichment/GO.${snp_set}.annotation

    # make significance tables with headers
    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/sig.out -d "," | \
         tr -d " \t"    | sed 's/FALSE/0/g' | sed 's/TRUE/1/g' \
         > ~/urchin_af/analysis/go_enrichment/sig_binary.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/logp.out -d "," | \
        tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/logp.${snp_set}.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/pval.out -d "," | \
         tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/pval.${snp_set}.table

    # change for topGO format, need to be comma delimted
    cat ~/urchin_af/analysis/go_enrichment/GO.${snp_set}.annotation | sed 's/;/,/g' | \
        grep -v 'unknown' > ~/urchin_af/analysis/go_enrichment/topGO.${snp_set}.annotation

done
