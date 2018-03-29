
# all variants
dat <- read.csv("~/urchin_af/analysis/cmh.master.out", header=TRUE, stringsAsFactor=FALSE, sep= "\t")
dat$logP_75 <- -log10(dat$qval_pH75)
dat$logP_80 <- -log10(dat$qval_pH80)

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
    a <- dat[which(dat$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_master_GO.out", out, col.names=TRUE,
    row.names=FALSE, quote=FALSE,sep="\t")

# pull out unique genes. ie, where multiple snps are in gene, only use that gene, rather than each variant
gene <- unique(dat$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(dat)
for(i in 1:length(gene)){
    a <- dat[which(dat$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_master_GO.out", out, col.names=TRUE,
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
    out[i,] <- a[which.min(as.numeric(a$qval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_genic_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# genic
gene <- unique(genic$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(genic)
for(i in 1:length(gene)){
    a <- genic[which(genic$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_genic_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



# exon
gene <- unique(exon$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(exon)
for(i in 1:length(gene)){
    a <- exon[which(exon$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_exon_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

gene <- unique(exon$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(exon)
for(i in 1:length(gene)){
    a <- exon[which(exon$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_exon_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


# intron
gene <- unique(intron$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(intron)
for(i in 1:length(gene)){
    a <- intron[which(intron$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_intron_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

gene <- unique(intron$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(intron)
for(i in 1:length(gene)){
    a <- intron[which(intron$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_intron_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# intergenic

gene <- unique(intergenic$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(intergenic)
for(i in 1:length(gene)){
    a <- intergenic[which(intergenic$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_intergenic_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

gene <- unique(intergenic$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(intergenic)
for(i in 1:length(gene)){
    a <- intergenic[which(intergenic$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_intergenic_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# non_coding
gene <- unique(non_coding$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(non_coding)
for(i in 1:length(gene)){
    a <- non_coding[which(non_coding$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH75)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH75_non_coding_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# non_coding
gene <- unique(non_coding$SPU_1)
out <- as.data.frame(matrix(nrow=length(gene), ncol=ncol(dat)))
colnames(out) <- colnames(non_coding)
for(i in 1:length(gene)){
    a <- non_coding[which(non_coding$SPU_1 == gene[i]),]
    out[i,] <- a[which.min(as.numeric(a$qval_pH80)),]
}

write.table(file="~/urchin_af/analysis/go_enrichment/cmh.pH80_non_coding_GO.out", out, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

