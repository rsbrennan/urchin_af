### Fig S4

# permutation results
# checking for bias in CMH stat

library(dplyr)
library(broom)
library(tidyr)
library(scales)
library(gridBase)
library(grid)
library(ggplot2)

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", stringsAsFactors=FALSE, header=TRUE)
bs <- read.table("~/urchin_af/analysis/permutation_af.txt", stringsAsFactors=FALSE, header=TRUE)
cut_off <- (0.05/9828)

mydata$D1_8_af <- (sapply(mydata$D1_8_af,function(x)  
          ifelse(x > 0.5, (1-x), x)))

# and fold all the permutations

(sapply(mydata$D1_8_af,function(x)  
          ifelse(x > 0.5, (1-x), x)))
bs[] <- lapply(bs, function(x)  
          ifelse(x > 0.5, (1-x), x))


# need to pull out only selected alleles
snp.sel_75 <- mydata$D1_8_af[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off)]
snp.sel_80 <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off)]
snp.sel_both <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off)]

# compare stats:


ph8 <- c()
ph7 <- c()
phboth <- c()


for (i in 1: ncol(bs)){
   ph8[i] <-  ks.test(snp.sel_80, bs[,i])$p.value
   ph7[i] <-  ks.test(snp.sel_75,bs[,i])$p.value
   phboth[i] <-  ks.test(snp.sel_both,bs[,i])$p.value
    #wilcox.test(snp.sel_80, bs[,i], alternative="l", correct=TRUE)$p.value,
    #wilcox.test(snp.sel_75, bs[,i], alternative="l", correct=TRUE)$p.value,
    #wilcox.test(snp.sel_both,bs[,i], alternative="l", correct=TRUE)$p.value
}

length(which(ph8 < 0.05))/500
length(which(ph7 < 0.05))/500
length(which(phboth < 0.05))/500

######
##
## Fig S4
##
######

# read in data from permutation in 06_balancing_sel_permutation.R

## plot CI

dat.bs <- data.frame(CHROM = mydata$CHROM, POS = mydata$POS, af_out = mydata$D1_8_af)

out <- list()
for (i in 1: ncol(bs)){

    out[[i]] <- data.frame(x=density(bs[,i], bw = 0.05)$x,
        y= density(bs[,i], bw = 0.05)$y,
        bin = cut(density(bs[,i], bw = 0.05)$x, breaks=seq(from=-0.3, to=1.3, by=0.01)))
}

out.new <- do.call(rbind, out)

densities.qtiles <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles$x <- seq(from=-0.15, to=0.65, by=0.01)

avg_perm <- c()
for (i in 1: ncol(bs)){
    avg_perm <- c(avg_perm, bs[,i])
}

png("~/urchin_af/figures/Fig_S4_cmhperm.png", height=100, width=100, units="mm", res=300)

par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(density(0:1), ylim=c(0,4),xlim=c(0,0.5), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Starting allele frequency", line=1.5, cex.lab=1)
title(ylab="Density", line=1.5, cex.lab=1)

polygon(x=c(densities.qtiles$x,rev(densities.qtiles$x)),
    y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
    col=alpha("black", alpha=0.1),border=NA)

lines(densities.qtiles$x,densities.qtiles$q50, col="black", lwd=3 )
#lines(density(mydata$D1_8_af), col=alpha("black", 1), lwd=3)
lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80,bw=0.05), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both,bw=0.05), col=alpha("darkorchid2", 1), lwd=3)

legend("topright", c("permuted", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=1,
    cex=0.8, lwd=2)

dev.off()

