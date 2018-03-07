# look for overlap between my analysis and pespeni and evans

dat <- read.csv("~/urchin_af/analysis/cmh.gene_names.out", header=TRUE, stringsAsFactors=FALSE, sep="\t")

dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_pval < 0.01)] <- FALSE

evans <- read.table("~/urchin_af/analysis/evans_list", header=FALSE)
pesp_icb <- read.table("~/urchin_af/analysis/pespeni_2013_icb_selected.txt", header=FALSE)
pesp_pnas <- read.table("~/urchin_af/analysis/pespeni_pnas.txt", header=FALSE)

dat_spu <- unique(unlist(strsplit(as.character(unique(dat$SPU_1)), "-")))

# identify overlap between gene lists.

pnas <- intersect(dat_spu, pesp_pnas$V1)
length(pnas)
# 18 of the 30 pnas genes are in my data

pnas_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(pnas)){
    out.tmp <- dat[grep(pnas[i], dat$SPU_1),]
    pnas_out <- rbind(pnas_out, out.tmp)
}

sig_pnas <- pnas_out[which(pnas_out$sig == TRUE),]
nrow(sig_pnas)
# 4 sig overlaps, 1 gene SPU_020412

# evans data
evans_spu <- intersect(dat_spu, evans$V1)
length(evans_spu)
#[1] 383
evans_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(evans_spu)){
    out.tmp <- dat[grep(evans_spu[i], dat$SPU_1),]
    evans_out <- rbind(evans_out, out.tmp)
}

sig_evans <- evans_out[which(evans_out$sig == TRUE),]
nrow(sig_evans)
# 57 of these significant evans genes are significant in my data
# 41 genes


### icb data
icb_spu <- intersect(dat_spu, pesp_icb$V1)
length(icb_spu)
#121 overlap

icb_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(icb_spu)){
    out.tmp <- dat[grep(icb_spu[i], dat$SPU_1),]
    icb_out <- rbind(icb_out, out.tmp)
}

sig_icb <- icb_out[which(icb_out$sig ==TRUE),]
nrow(sig_icb)
# 11 overlapping genes
unique(sig_icb$SPU_1)
# 10 unique genes

write.table(file="~/urchin_af/analysis/icb_overlap.txt"  ,sig_icb, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
write.table(file="~/urchin_af/analysis/pnas_overlap.txt" ,sig_pnas, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
write.table(file="~/urchin_af/analysis/evans_overlap.txt",sig_evans, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


##############################
##### test for sig overlap
##############################

pesp_all <- read.table("~/urchin_af/variants/Pespeni_LarvalSelection_spu.txt",
    header=TRUE, stringsAsFactors=FALSE)
# pull out significant spu from cmh

dat.sig <- dat[which(dat$sig == TRUE),]
sig_spu <- unique(unlist(strsplit(as.character(unique(dat.sig$SPU_1)), "-")))

# all possible genes included
pnas <- unique(pesp_all$X.CHROM)
universe <- length(unique(c(dat_spu, pnas)))
sig_pnas <- as.character(pesp_pnas$V1)

mat <- matrix(
   c( length(intersect(sig_spu, sig_pnas)),
      length(setdiff(sig_pnas, sig_spu)),
      length(setdiff(sig_spu, sig_pnas)),
      universe - length(union(sig_spu, sig_pnas))),
   nrow=2
   )

fr <- fisher.test(mat, alternative="greater")
fr

# all possible genes included
icb_sig <- as.character(pesp_icb$V1)

pnas <- unique(pesp_all$X.CHROM)
universe <- length( unique(c(dat_spu, pnas)))

mat <- matrix(
   c( length(intersect(sig_spu, icb_sig)),
      length(setdiff(icb_sig, sig_spu)),
      length(setdiff(sig_spu, icb_sig)),
      universe - length(union(sig_spu, icb_sig))),
   nrow=2
   )

fr <- fisher.test(mat, alternative="greater")
fr

# all possible genes included
evans_sig <- as.character(evans$V1)

pnas <- unique(pesp_all$X.CHROM)
universe <- length( unique(c(dat_spu, pnas)))

mat <- matrix(
   c( length(intersect(sig_spu, evans_sig)),
      length(setdiff(evans_sig, sig_spu)),
      length(setdiff(sig_spu, evans_sig)),
      universe - length(union(sig_spu, evans_sig))),
   nrow=2
   )

fr <- fisher.test(mat, alternative="greater")
fr
