
# all variants
dat <- read.csv("~/urchin_af/analysis/cmh.master.out", header=TRUE, stringsAsFactor=FALSE, sep= "\t")
dat$logP <- -log10(dat$PVAL)

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
    a <- dat[which(dat$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.master_GO.out", out, col.names=TRUE, 
    row.names=FALSE, quote=FALSE,sep="\t")

##############
# subset snps to different groups

exon <- dat[which(dat$class == "synonymous" | dat$class == "non-synonymous"),]
intron <- dat[which(dat$class == "intron"),]
intergenic <- dat[which(dat$class == "intergenic"),]
synonymous <- dat[which(dat$class == "synonymous"),]
non_syn <- dat[which(dat$class == "non-synonymous"),]
non_coding <- dat[which(dat$class == "intergenic" | dat$class == "intron"),]
genic <- dat[which(dat$class == "synonymous" | dat$class == "non-synonymous" | dat$class == "intron"),]

# parse down to remove multiple snps in genes.

# genic
gene <- unique(genic$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(genic)
for(i in 1:length(gene)){
    a <- genic[which(genic$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.genic_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# exon
gene <- unique(exon$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(exon)
for(i in 1:length(gene)){
    a <- exon[which(exon$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.exon_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# intron
gene <- unique(intron$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(intron)
for(i in 1:length(gene)){
    a <- intron[which(intron$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.intron_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# intergenic
gene <- unique(intergenic$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(intergenic)
for(i in 1:length(gene)){
    a <- intergenic[which(intergenic$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.intergenic_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# synonymous
gene <- unique(synonymous$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(synonymous)
for(i in 1:length(gene)){
    a <- synonymous[which(synonymous$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.synonymous_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# non_syn
gene <- unique(non_syn$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(non_syn)
for(i in 1:length(gene)){
    a <- non_syn[which(non_syn$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.non_syn_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# non_coding
gene <- unique(non_coding$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(non_coding)
for(i in 1:length(gene)){
    a <- non_coding[which(non_coding$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$PVAL)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.non_coding_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
