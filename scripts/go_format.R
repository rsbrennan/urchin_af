
# all variants
dat <- read.table("~/urchin_af/analysis/cmh.master.out", header=TRUE)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$WHL)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
    a <- dat[which(dat$WHL == gene[i]),]
    out[i,] <- a[which.min(a$PVAL),]
}

out$logP <- -log10(out$PVAL)

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.masterGO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# coding
# read in exon and nc list
exon <- read.table("~/urchin_af/analysis/variants_exon.bed", header=FALSE)
nonc <- read.table("~/urchin_af/analysis/variants_noncoding.bed", header=FALSE)

# convert to snps
exon$SNP <- paste(exon$V1, exon$V3, sep=":")
nonc$SNP <- paste(nonc$V1, nonc$V3, sep=":")
exonSel$SNP <- paste(exonSel$V1, exonSel$V3, sep=":")
noncSel$SNP <- paste(noncSel$V1, noncSel$V3, sep=":")

#from dat, pull out coding from entire list

exon_dat<- dat[dat$SNP %in% exon$SNP,] 
nonc_dat<- dat[dat$SNP %in% nonc$SNP,] 

# parse down to remove multiple snps in genes. EXONS first
i <- sapply(exon_dat, is.factor)
exon_dat[i] <- lapply(exon_dat[i], as.character)
gene <- unique(exon_dat$WHL)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(exon_dat)))
colnames(out) <- colnames(exon_dat)
for(i in 1:length(gene)){
    a <- exon_dat[which(exon_dat$WHL == gene[i]),]
    out[i,] <- a[which.min(a$PVAL),]
}

out$logP <- -log10(out$PVAL)

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.codingGO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# parse down to remove multiple snps in genes. NONCODING next
i <- sapply(nonc_dat, is.factor)
nonc_dat[i] <- lapply(nonc_dat[i], as.character)
gene <- unique(nonc_dat$WHL)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(nonc_dat)))
colnames(out) <- colnames(nonc_dat)
for(i in 1:length(gene)){
    a <- nonc_dat[which(nonc_dat$WHL == gene[i]),]
    out[i,] <- a[which.min(a$PVAL),]
}

out$logP <- -log10(out$PVAL)

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.NONcodingGO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
