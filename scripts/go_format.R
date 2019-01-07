# go_format.R

# all variants
dat <- read.csv("~/urchin_af/analysis/cmh.master.out", header=TRUE, stringsAsFactor=FALSE, sep= "\t")
dat$logP_75 <- -log10(dat$pval_pH75)
dat$logP_80 <- -log10(dat$pval_pH80)
dat$sig_pH75 <- gsub(" TRUE", "TRUE",dat$sig_pH75)
dat$sig_pH80 <- gsub(" TRUE", "TRUE",dat$sig_pH80)
# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
    a <- dat[which(dat$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$pval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_master_GO.out", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep="\t")

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
    a <- dat[which(dat$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$pval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_master_GO.out", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep="\t")


#### overlap between the two ph
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
out$overlap_sig <- FALSE

for(i in 1:length(gene)){
    a <- dat[which(dat$SPU_1 == gene[i]),]
    if(sum(a$sig_pH75 == TRUE) >= 1 & sum(a$sig_pH80 == TRUE) >= 1){
        out[i,] <- a[which.min(as.numeric(a$pval_pH80)),]
        out$overlap_sig[i] <- TRUE
    }
    if(sum(a$sig_pH75 == TRUE) >= 1 & sum(a$sig_pH80 == TRUE) == 0){
        out[i,] <- a[which.min(as.numeric(a$pval_pH75)),]
        out$overlap_sig[i] <- FALSE
    }
    if(sum(a$sig_pH75 == TRUE) == 0 & sum(a$sig_pH80 == TRUE) >= 1){
        out[i,] <- a[which.min(as.numeric(a$pval_pH80)),]
        out$overlap_sig[i] <- FALSE
    }
    if(sum(a$sig_pH75 == TRUE) == 0 & sum(a$sig_pH80 == TRUE) == 0){
        out[i,] <- a[which.min(as.numeric(a$pval_pH75)),]
        out$overlap_sig[i] <- FALSE
    }
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.overlap_master_GO.out", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep="\t")

