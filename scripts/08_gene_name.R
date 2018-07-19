### the gene descriptions are terrible in some cases. merge with old annotation to try to improve them

dat <- read.csv("~/urchin_af/analysis/cmh.master.out", header=TRUE, sep="\t")
g.nm <- read.csv("~/urchin_af/analysis/cmh.master.whl_old.txt", header=TRUE, sep="\t")

dat$old_name <- g.nm$gene_name[match(dat$SNP_1, g.nm$SNP)]

write.table(file = "~/urchin_af/analysis/cmh.gene_names.out", dat,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
