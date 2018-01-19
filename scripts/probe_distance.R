library(scales)
dat <- read.table("~/urchin_af/analysis/snp_probe_dist.txt", header=FALSE)

af.dat <- read.table("~/urchin_af/data/allele.freq_all.txt", header=TRUE)

dat$SNP <- paste(dat$V1, dat$V3, sep=":")
af.dat$SNP <- paste(af.dat$CHROM, af.dat$POS, sep=":")

mydat <- merge(dat, af.dat, by="SNP")

cov.h <- mydat[which(mydat$V14 > 2000),]
cov.l <- mydat[which(mydat$V14 <= 2000),]

dep.h <-  apply(cov.h[,grep("_DPtotal", colnames(cov.h))], 1,mean)   
dep.l <-  apply(cov.l[,grep("_DPtotal", colnames(cov.l))], 1,mean)   

png("~/urchin_af/figures/depth_probe_dist.png", res=300, height=7, width=7, units="in")
hist(dep.l, breaks=40, col='grey', freq=FALSE, ylim=c(0,0.017), xlab="average depth of coverage", main="On- vs. Off-target avg depth")
hist(dep.h, breaks=30, col=alpha('red', 0.5), add=T, freq=FALSE)
legend(x=200, y= 0.014, c("On-target", "Off-target"), pch=15, col=c("grey","red"), cex=1.4, pt.cex=2)

dev.off()

###
### remove loci > 2kb from any target
###

# order both by SNP
dat.new <- dat[order(dat$SNP),]
af.dat.new <- af.dat[order(af.dat$SNP),]

low.keep <- which(dat.new$V14 <= 2000)

out <- af.dat.new[low.keep,]

write.table(file="~/urchin_af/data/allele.freq.txt",out,row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# write file to filter vcf
write.table(file="~/urchin_af/variants/keep.ontarget.txt", 
    as.data.frame(cbind(as.character(out$CHROM), out$POS)), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
