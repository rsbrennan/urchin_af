# look for overlap between pH 7.5 and 8.0- at the gene level.

dat <- read.csv("~/urchin_af/analysis/cmh.master_goodnm.out", header=TRUE, stringsAsFactors=FALSE, sep="\t")

dat$sig_pH75[which(dat$qval_pH75 < 0.001)] <- TRUE
dat$sig_pH80[which(dat$qval_pH80 < 0.001)] <- TRUE

sig_pH75 <- dat[which(dat$qval_pH75 < 0.001),]
sig_pH80 <- dat[which(dat$qval_pH80 < 0.001),]
sig_overlap <- dat[which(dat$qval_pH80 < 0.001 & dat$qval_pH75 < 0.001),]


pH75_spu <- unique(unlist(strsplit(as.character(unique(sig_pH75$SPU_1)), "-")))
pH80_spu <- unique(unlist(strsplit(as.character(unique(sig_pH80$SPU_1)), "-")))
overlap_spu <- unique(unlist(strsplit(as.character(unique(sig_overlap$SPU_1)), "-")))

pH75_spu <- pH75_spu[grep ("CHR_", pH75_spu, invert=TRUE)]
pH80_spu <- pH80_spu[grep ("CHR_", pH80_spu, invert=TRUE)]

# identify overlap between gene lists.

overlap <- intersect(pH75_spu, pH80_spu)
length(overlap)
# 207

overlap_out <- as.data.frame(matrix(ncol=ncol(dat),nrow=0))
for(i in 1:length(overlap)){
    out.tmp <- dat[grep(overlap[i], dat$SPU_1),]
    overlap_out <- rbind(overlap_out, out.tmp)
}

sig_overlap <- overlap_out[which(overlap_out$sig_pH80 == TRUE | overlap_out$sig_pH75 == TRUE ),]


write.table(file="~/urchin_af/analysis/gene_overlap.txt", sig_overlap,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
