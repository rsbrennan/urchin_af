
for snp_set in masterGO NONcodingGO codingGO; do

    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}.out | cut -f 1 | tail -n +2 > ~/urchin_af/analysis/go_enrichment/snp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}.out | cut -f 1  > ~/urchin_af/analysis/go_enrichment/snp_head.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}.out | cut -f 9 | tail -n +2 > ~/urchin_af/analysis/go_enrichment/GO.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}.out | cut -f 5  > ~/urchin_af/analysis/go_enrichment/sig.out

    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}.out | cut -f 10  > ~/urchin_af/analysis/go_enrichment/logp.out
    cat ~/urchin_af/analysis/go_enrichment/cmh.${snp_set}.out | cut -f 4  > ~/urchin_af/analysis/go_enrichment/pval.out

    # make go annotation
    paste ~/urchin_af/analysis/go_enrichment/snp.out  ~/urchin_af/analysis/go_enrichment/GO.out >  ~/urchin_af/analysis/go_enrichment/GO.${snp_set}.annotation

    # make significance tables
    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/sig.out -d "," | tr -d " \t"    | sed 's/FALSE/0/g' | sed 's/TRUE/1/g' > ~/urchin_af/analysis/go_enrichment/sig_binary.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/logp.out -d "," | tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/logp.table

    paste ~/urchin_af/analysis/go_enrichment/snp_head.out ~/urchin_af/analysis/go_enrichment/pval.out -d "," | tr -d " \t  "  > ~/urchin_af/analysis/go_enrichment/pval.table

    # change for topGO format, need to be comma delimted

    cat ~/urchin_af/analysis/go_enrichment/GO.${snp_set}.annotation | sed 's/;/,/g' | grep -v 'unknown' >~/urchin_af/analysis/go_enrichment/topGO.${snp_set}.annotation

done

