# Purple Sea Urchin Sequence Pool-seq Capture Selection Experiment

This repository holds scripts for the analysis of a purple sea urchin, Strongylocentrotus purpuratus, selection experiment.

Original experiment was completed on ?? where 25 adult urchins were spawned and resulting embryos were exposed to low pH or control conditions. Pools of larvae sequenced to > 50x coverage to identify adaptive changes in allele frequency in response to low pH selection.

## data availability

Sequence data: link coming soon.

GO annotations were retrived from ensembl on 2/18/18 and exported though biomart. This correponds to assembly Spur_3.1. http://metazoa.ensembl.org/Strongylocentrotus_purpuratus/Info/Index


## programs used

- SnpEff Version 4.3 build 2017-11-24
- R version 3.4.3
    + pcadapt 3.0.4
- VCFtools version 0.1.15 
- LDx
- topGO version 2.22.0 
- Python 2.7.5

Program versions can also be found in the scripts and/or in the associated manuscript (currently in prep)

## Scripts

### Alignment, variant calling, filtering

Get alignment, trimming script from april. 

Call variants: `01.1_variant_call.sh`

#### filtering

- Initial filtering for depth, bi-allelic: `01.2_vcf_filter.sh`
- Pull out average depth. To assign max depth filter: `01.3_filter.R`
- Filter max mean depth of 372: `01.4_filter_depth.sh`
- Off target filtering
    - Find distance to closest probe: `01.5_probe_distance.sh`
    - Output those within 2kb: `01.6_probe_distance.R`
    - keep on target: `01.7_filter_final.sh`

### Analysis

- identify allele frequency changes: `02_AF_change.R`
- polarize by "adaptive" allele: `03_polarize_adaptive.R`
- assign overlap with gene models: `04_snpeff.sh`
- Look at allele frequencies of selected loci: `05_balancing_sel.R`
- Compare af distribution of selected loci vs 1% quantile of permuted samples. Basically to ask if we are pulling these out only because our cmh analysis is biased
    - `06_balancing_sel_permutation.R`
- Assign GO terms to SNPs: `07_snp_go_assign.py`
- GO enrichment: `09_go_enrich.R`
- LD decay: `11_ld.sh`
- Look for overlap of selected genes with old studies: `12_gene_comparison_old_studies.R`


### Figures

- Figure 1:
- Figure 2: `Fig_02_pca.R`
- Figure 3: `Fig_03_LD.R`
- Figure 4: `Fig_04_afchangecomp.R`
- Figure 5: `Fig_05.R`

Supplemental Figures:
- Figure S1: in `Fig_03_LD.R`
- Figure S2: in `Fig_03_LD.R`
- Figure S3: `Fig_S2_manhattan.R`
- Figure S4:
- Figure S5: 
- Figure S6: 




