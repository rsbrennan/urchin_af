# look for overlap between my analysis and pespeni and evans

dat <- read.table("~/urchin_af/analysis/cmh.master.out", header=TRUE)

evans <- read.table("~/urchin_af/analysis/evans_list", header=FALSE)
pesp_icb <- read.table("~/urchin_af/analysis/pespeni_2013_icb_selected.txt", header=FALSE)
pesp_pnas <- read.table("~/urchin_af/analysis/pespeni_pnas.txt", header=FALSE)


dat_spu <- unique(unlist(strsplit(as.character(unique(dat$SPU)), ";")))

# identify overlap between gene lists.

pnas <- intersect(dat_spu, pesp_pnas$V1)
length(pnas)
# 16 of the 30 pnas genes are in my data

pnas_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(pnas)){
    out.tmp <- dat[grep(pnas[i], dat$SPU),]
    pnas_out <- rbind(out, out.tmp)
}

sig_pnas <- pnas_out[which(pnas_out$sig ==TRUE),]
nrow(sig_pnas)
# only 1 significant SPU overlaps: translation_initiation_factor_2_gamma_subunit WHL22.669828 SPU_020412

# evans data
evans_spu <- intersect(dat_spu, evans$V1)

evans_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(evans_spu)){
    out.tmp <- dat[grep(evans_spu[i], dat$SPU),]
    evans_out <- rbind(evans_out, out.tmp)
}

sig_evans <- evans_out[which(evans_out$sig ==TRUE),]
nrow(sig_evans)
# 27 of these significant evans genes are significant in my data


### icb data
icb_spu <- intersect(dat_spu, pesp_icb$V1)
#120 overlap

icb_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(icb_spu)){
    out.tmp <- dat[grep(icb_spu[i], dat$SPU),]
    icb_out <- rbind(icb_out, out.tmp)
}

sig_icb <- icb_out[which(icb_out$sig ==TRUE),]
unique(sig_icb$gene_name)
# 6 overlapping genes


###### look at functional annotation of these snps

ann <- read.table("~/urchin_af/analysis/cmh.annotations.out", header=TRUE, sep=" ")

ann_icb <- merge(ann, sig_icb, by="SNP")
ann_evans <- merge(ann, sig_evans, by="SNP")
