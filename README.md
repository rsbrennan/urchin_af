# Purple Sea Urchin Sequence Pool-seq Capture Selection Experiment

This repository holds scripts for the analysis of a purple sea urchin, *Strongylocentrotus purpuratus*, selection experiment. Now published in Proceedings of the Royal Society B. [Link to paper here](https://doi.org/10.1098/rspb.2019.0943). 


If you use these resources, please cite the following publication:

Brennan, R.S., Garrett, A.D., Huber, K.E., Hargarten, H. and Pespeni, M.H., 2019. Rare genetic variation and balanced polymorphisms are important for survival in global change conditions. Proceedings of the Royal Society B, 286(1904), p.20190943.

Original experiment was run where 25 adult urchins were spawned and resulting embryos were exposed to moderately low (pH 8.0) or extremely low (pH 7.5) pH conditions. Pools of larvae sequenced to > 50x coverage to identify adaptive changes in allele frequency in response to different degrees of low pH selection.

## data availability

Sequence data: raw reads are available at NCBI SRA: BioProject ID: PRJNA479817

GO annotations were retrived from ensembl on 2/18/18 and exported though biomart. This correponds to assembly Spur_3.1. http://metazoa.ensembl.org/Strongylocentrotus_purpuratus/Info/Index


## programs used

- SnpEff Version 4.3 build 2017-11-24
- R version 3.4.3
    + pcadapt 3.0.4
- VCFtools version 0.1.15
- [LDx](https://sourceforge.net/p/ldx/wiki/Home/)
- topGO version 2.22.0 
- Python 2.7.5

Program versions can also be found in the scripts and/or in the associated manuscript (currently in prep)

## Scripts

### Alignment, variant calling, filtering

Alignment, trimming: `01.0_CleanAndMap.sh`  

Call variants: `01.1_variant_call.sh`

#### filtering

- Initial filtering for depth, bi-allelic: 
    - `01.2_vcf_filter.sh`
- Pull out average depth. To assign max depth filter: 
    - `01.3_filter.R`
- Filter max mean depth of 372: 
    - `01.4_filter_depth.sh`
- Off target filtering
    - Find distance to closest probe: 
        - `01.5_probe_distance.sh`
    - Output those within 2kb: 
        - `01.6_probe_distance.R`
    - keep on target: 
        - `01.7_filter_final.sh`

### Analysis

- identify allele frequency changes with CMH: 
    - `02_AF_change.R`
- polarize by "adaptive" allele: 
    - `03_polarize_adaptive.R`
- assign overlap with gene models: 
    - `04_snpeff.sh`
- Look at allele frequencies of selected loci: 
    - `05_balancing_sel.R`
- Compare af distribution of selected loci vs permuted samples. To ask if we are pulling these out only because our cmh analysis is biased
    - `06_balancing_sel_permutation.R`
- Assign GO terms to SNPs: 
    - `07_snp_go_assign.py`
    - some assigned names aren't very informative. Use transcriptome to update these or fill in missing
        - `whl_assign.sh`
        - `whl_assign.py`
        - these output ``cmh.master_goodnm.out`, which is the final output
- GO enrichment: 
    - Choose only 1 snp per gene: 
        - `go_format.R`
    - then format: 
        - `go_format.sh`
    - run enrichment: 
        - `09_go_enrich.R`
- LD decay: 
    - `11_ld.sh` (I know, I know, there's no #10)
- Look for overlap of selected genes with old studies: 
    - `12_gene_comparison_old_studies.R`
- PCA: 
    - `Fig_02_pca.R`
- Nucleotide diversity: 
    - first, sam to bam, then run popoolation. 
        - `pi_popoolation.sh`
    - Then calc with R: 
        - `pi_popoolation.R`
- Morphology
    - `morphology.R`


### Figures

- Figure 1: 
- Figure 2: `Fig_02.R`
- Figure 3: `Fig_03.R`
- Figure 4: `Fig_04.R`

Supplemental Figures:
- Figure S1: `Fig_S1.R`
- Figure S2: `Fig_S02_03_05.R`
- Figure S3: `Fig_S02_03_05.R`
- Figure S4: `Fig_S04.R` 
- Figure S5: `Fig_S02_03_05.R`
- Figure S6: `05_balancing_sel.R`




