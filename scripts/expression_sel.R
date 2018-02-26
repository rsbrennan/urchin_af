
library(scales)
# looking at gene expression and selection

adult <- read.csv("~/urchin_af/data/Pespeni_SpurpuratusAdultCommonGardenGeneExpressionData.txt",
    header=TRUE, sep="\t")

embryo <- read.csv("~/urchin_af/data/CombinedCounts_NotNormalized.txt",
    header=TRUE, sep="\t")

a_mean <- data.frame(ID=adult[,1], Means=rowMeans(adult[,-1]))
e_mean <- data.frame(ID=embryo[,1], Means=rowMeans(embryo[,-1]))

dat <- merge(a_mean, e_mean, by="ID")
colnames(dat) <- c("spu", "adult", "embryo")

# read in cmh data

cmh <- read.csv("~/urchin_af/analysis/cmh.master.out", header=TRUE, sep='\t')

new <- merge(dat, cmh, by.x="spu", by.y="SPU_1")

new.sig <- new[which(new$sig == TRUE),]

which(new$embryo == Inf)
dat$embryo[3978] <- mean(dat$embryo)
dat$embryo[901] <- mean(dat$embryo)
dat$embryo[which(dat$embryo ==0)] <- mean(dat$embryo)

mean(log10(new.sig$adult))
mean(log10(new.sig$embryo))

median(log10(new.sig$adult))
median(log10(new.sig$embryo))

png("~/urchin_af/figures/expression_selection_norm.png", 
    height=7, width=7, res=300, units="in")
plot(y=(log10(dat$adult)/mean(log10(dat$adult))), 
    x=(log10(dat$embryo)/mean(log10(dat$embryo))),
    col=alpha("black", 0.1), pch=19,
    xlab="log10 mean embryo expression",
    ylab="log10 mean adult expression",
    xlim=c(0,3), ylim=c(0,3))
points(y=log10(new.sig$adult)/mean(log10(dat$adult)), 
    x=log10(new.sig$embryo)/mean(log10(dat$embryo)),col="red", pch=19 )

dev.off()


mean(log10(new.sig$adult)/mean(log10(dat$adult)))
mean(log10(new.sig$embryo)/mean(log10(dat$adult)))
 
median(log10(new.sig$adult))
median(log10(new.sig$embryo))
