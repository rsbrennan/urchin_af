mydat <- read.table("~/urchin_af/analysis/snp_probe_dist.txt", header=FALSE)

###
### remove loci > 2kb from any target
###

low.keep <- which(mydat$V14 <= 2000)

out <- mydat[low.keep,]

# write file to filter vcf
write.table(file="~/urchin_af/variants/keep.ontarget.txt",
    as.data.frame(cbind(as.character(out$V1), out$V3)), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
