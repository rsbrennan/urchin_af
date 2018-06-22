library(stringr)
library(ggplot2)
library(gridExtra)
library(MASS)

dat <- read.csv("~/urchin_af/analysis/cmh.master.out", 
    stringsAsFactors=FALSE, header=TRUE, sep="\t")

afdata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE, stringsAsFactors=FALSE)

afdata$SNP <- paste(afdata$CHROM, afdata$POS, sep=":")

new <- merge(dat, afdata, by.x="SNP_1", by.y="SNP")

sig_pH75 <- new[which(new$pH7_selection_qval < 0.001),]
sig_pH80 <- new[which(new$pH8_selection_qval < 0.001),]

ph7_exon <- sig_pH75[which(sig_pH75$class == "synonymous" | 
                    sig_pH75$class == "non-synonymous"),]

ph7_nc <- sig_pH75[which(sig_pH75$class == "intergenic" | 
                    sig_pH75$class == "intron"),]

ph7_exon$af_change <- abs(ph7_exon$D7_7_mean - ph7_exon$D1_8_mean)
ph7_nc$af_change <- abs(ph7_nc$D7_7_mean - ph7_nc$D1_8_mean)



ph8_exon <- sig_pH80[which(sig_pH80$class == "synonymous" | 
                    sig_pH80$class == "non-synonymous"),]

ph8_nc <- sig_pH80[which(sig_pH80$class == "intergenic" | 
                    sig_pH80$class == "intron"),]

ph8_exon$af_change <- abs(ph8_exon$D7_8_mean - ph8_exon$D1_8_mean)
ph8_nc$af_change <- abs(ph8_nc$D7_8_mean - ph8_nc$D1_8_mean)


png("~/urchin_af/figures/coding_vs_non.png", height=100, width=200, units="mm", res=300)
par(mfrow = c(1, 2), mar=c(4, 4, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(density(ph8_nc$af_change, bw=0.01), col="red", lwd=2,
    main="pH 8.0",xlab="mean AF change",
    ylim=c(-0.5,10), xlim=c(0, 0.4))
lines(density(ph8_exon$af_change, bw=0.01), col="black", lwd=2)

rug(ph8_exon$af_change)
rug(ph8_nc$af_change, col="red", line=-1)

plot(density(ph7_nc$af_change, bw=0.01), col="red", lwd=2,
    main="pH 7.5", xlab="mean AF change", ylab="",
    ylim=c(-0.5,10), xlim=c(0, 0.4))
lines(density(ph7_exon$af_change, bw=0.01), col="black", lwd=2)

rug(ph7_exon$af_change)
rug(ph7_nc$af_change, col="red", line=-1)

dev.off()

