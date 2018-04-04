

# looking at gene expression and selection

adult <- read.csv("~/urchin_af/data/UniqueCounts12Individuals_e1o2k8.txt",
    header=TRUE, sep="\t")

embryo <- read.csv("~/urchin_af/data/CombinedCounts_NotNormalized.txt",
    header=TRUE, sep="\t")

colnames(adult) <- c("ID",paste("adult",str_split_fixed(colnames(adult)[2:ncol(adult)], "_", 2)[,2], sep="_"))

# merge data
dat <- merge(adult, embryo, by.x="ID", by.y="GeneName")
rownames(dat) <- dat$ID
dat <- dat[2:ncol(dat)]

# need a file that specifies treatment, etc for each indiv

condition <- c(rep("0385", 12), str_split_fixed(colnames(dat)[13:ncol(dat)], "_", 5)[,4])

devstage <- c(rep("adult", 12), rep("larvae", 28))
#devstage <- c(rep("adult", 12), str_split_fixed(colnames(dat)[13:ncol(dat)], "_", 5)[,4])

coldata <- data.frame(condition=condition, devstage=devstage)
rownames(coldata) <- colnames(dat)

# convert to deseq object
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = coldata,
                              design = ~ devstage)

dds <- DESeq(dds)
res <- results(dds, alpha=0.01)
res

summary(res)

# out put specific comparisons
#res.A.B <- results(dds1, contrast=c("devstage","A","B"))


################  Data visualization

resLFC <- lfcShrink(dds, coef="devstage_larvae_vs_adult")

plotMA(res, main="DESeq2", ylim=c(-12,10))
plotMA(resLFC, main="DESeq2", ylim=c(-12,10))
#abline(h=c(-1,1), col="blue", lwd =2)

rld <- rlog(dds, blind=FALSE)
vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd,intgroup=c("devstage"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))


ggplot(data, aes(PC1, PC2, color=devstage, shape=devstage)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")


resSig <- subset(resOrdered, padj < 0.1)
resSig

# take mean of each treatment

baseMeanPerLvl <- sapply(levels(dds$devstage), function(lvl) 
    rowMeans( counts(dds,normalized=TRUE)[,dds$devstage == lvl] ) )
baseMeanPerLvl <- as.data.frame(baseMeanPerLvl)
baseMeanPerLvl$foldch <- res$log2FoldChange
baseMeanPerLvl$SPU <- rownames(baseMeanPerLvl)

# read in cmh data
cmh <- read.csv("~/urchin_af/analysis/cmh.master.out", header=TRUE, sep='\t')

new <- merge(baseMeanPerLvl, cmh, by.x="SPU", by.y="SPU_1")

new.sig <- new[which(new$PVAL < 0.0001),]

png("~/urchin_af/figures/expression_selection.png",
    height=7, width=7, res=300, units="in")
plot(y=log2(new$larvae),
    x=log2(new$adult),
    col=alpha("black", 0.1), pch=19,
    xlab="mean adult expression",
    ylab="mean embryo expression")
points(y=log2(new.sig$larvae),
    x=log2(new.sig$adult),col=alpha("red", 0.7), pch=19 )

dev.off()

png("~/urchin_af/figures/expression_selection_hist.png",
    height=7, width=10, res=300, units="in")

par(mfrow=c(1,2)) 

hist(log2(new$larvae), col=alpha("black", 0.2), freq=F, main="Larvae")
hist(log2(new.sig$larvae), col=alpha("red", 0.3), breaks=30, freq=F, add=T)

hist(log2(new$adult), col=alpha("black", 0.2), freq=F, main="adult")
hist(log2(new.sig$adult), col=alpha("red", 0.3), breaks=20, freq=F, add=T)

dev.off()

mean(log10(new.sig$adult)/mean(log10(dat$adult)))
mean(log10(new.sig$embryo)/mean(log10(dat$adult)))

median(log10(new.sig$adult))
median(log10(new.sig$embryo))
