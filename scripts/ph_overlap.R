# look for overlap between pH 7.5 and 8.0- at the gene level.

dat <- read.csv("~/urchin_af/analysis/cmh.annotations.out", header=TRUE, stringsAsFactors=FALSE, sep="\t")

length(unique(unlist(strsplit(as.character(unique(dat$Gene_ID)), "-"))))


sig_pH75 <- dat[which(dat$sig_pH7 == TRUE),]
nrow(sig_pH75)
#[1] 787
sig_pH80 <- dat[which(dat$sig_pH8 == TRUE),]
nrow(sig_pH80)
#[1] 572

sig_overlap <- dat[which(dat$sig_pH8 == TRUE & dat$sig_pH7 == TRUE),]
nrow(sig_overlap)
# [1] 128

pH75_spu <- unique(unlist(strsplit(as.character(unique(sig_pH75$Gene_ID)), "-")))
pH80_spu <- unique(unlist(strsplit(as.character(unique(sig_pH80$Gene_ID)), "-")))
overlap_spu <- unique(unlist(strsplit(as.character(unique(sig_overlap$Gene_ID)), "-")))

pH75_spu <- pH75_spu[grep ("CHR_", pH75_spu, invert=TRUE)]
pH80_spu <- pH80_spu[grep ("CHR_", pH80_spu, invert=TRUE)]
overlap_spu <- overlap_spu[grep ("CHR_", overlap_spu, invert=TRUE)]

length(pH75_spu)
#[1] 606
length(pH80_spu)
#[1] 470
length(overlap_spu)
#100

# identify overlap between gene lists.

overlap <- intersect(pH75_spu, pH80_spu)
length(overlap)
# 138

# 138/938

overlap_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(overlap)){
    out.tmp <- dat[grep(overlap[i], dat$Gene_ID),]
    overlap_out <- rbind(overlap_out, out.tmp)
}

sig_overlap <- overlap_out[which(overlap_out$sig_pH8 == TRUE | overlap_out$sig_pH7 == TRUE ),]


write.table(file="~/urchin_af/analysis/gene_overlap.txt", sig_overlap,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# count up total sig snps:
length(unique(c(sig_pH75$SNP, sig_pH80$SNP)))
#[1] 1231
1231/75368
#[1] 0.01633319

# total sig genes:
length(unique(c(pH75_spu, pH80_spu)))
#938
938/9828
#[1] 0.0954416

