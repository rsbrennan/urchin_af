### Fig S3

# permutation results
# checking for bias in CMH stat

library(dplyr)
library(broom)
library(tidyr)
library(scales)



mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", stringsAsFactors=FALSE, header=TRUE)
perm_dat <- read.table("~/urchin_af/analysis/permutation_af.txt", stringsAsFactors=FALSE, header=TRUE)
cut_off <- 0.001
# need to pull out only selected alleles
snp.sel_75 <- mydata$af_out[which(mydata$pH7_selection_qval < cut_off & mydata$pH8_selection_qval >= cut_off)]
snp.sel_80 <- mydata$af_out[which(mydata$pH8_selection_qval < cut_off & mydata$pH7_selection_qval >= cut_off)]
snp.sel_both <- mydata$af_out[which(mydata$pH8_selection_qval < cut_off & mydata$pH7_selection_qval < cut_off)]

png("~/urchin_af/figures/Fig_S3_afchangecomp.png", height=100, width=200, units="mm", res=300)
par(mfrow = c(1, 2), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

######
##
## maf plot
##
######

# read in data from permutation in 06_balancing_sel_permutation.R


plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.6, tcl=-0.2)
title(ylab="Density", line=1.5, cex.lab=0.7)
title(xlab="Starting allele frequency", line=1.5, cex.lab=0.7)

avg_perm <- c()
for (i in 1: ncol(perm_dat)){
    lines(density(unlist(perm_dat[,i]),na.rm=TRUE, bw=0.05), col=alpha("black", 0.08), lwd=1)
    avg_perm <- c(avg_perm, unlist(perm_dat[,i]))
}

lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80, bw=0.05), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both, bw=0.05), col=alpha("darkorchid2", 1), lwd=3)
#lines(density(avg_perm), col=alpha("darkgoldenrod3", 1), lwd=3, lty=2)

#abline(v=mean(avg_perm), lty=2, col= "black", lwd=3)
#abline(v=mean(snp.sel_75), lty=2, col= "firebrick3", lwd=3)
#abline(v=mean(snp.sel_80), lty=2, col= "royalblue3", lwd=3)
#abline(v=mean(snp.sel_both), lty=2, col= "darkorchid2", lwd=3)

legend("topright", c("permuted", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=1,
    cex=0.7, lwd=2)


mtext(text="A",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)


## plot CI

dat.bs <- data.frame(CHROM = mydata$CHROM, POS = mydata$POS, af_out = mydata$af_out)

bs <- read.table("~/urchin_af/analysis/permutation_af.txt")

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

densities.qtiles$x <- seq(from=-0.155, to=1.135, by=0.01)


avg_perm <- c()
for (i in 1: ncol(bs)){
    avg_perm <- c(avg_perm, bs[,i])
}


plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.6, tcl=-0.2)
title(xlab="Starting allele frequency", line=1.5, cex.lab=0.7)

polygon(x=c(densities.qtiles$x,rev(densities.qtiles$x)),
    y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
    col=alpha("black", alpha=0.1),border=NA)

lines(densities.qtiles$x,densities.qtiles$q50, col="black", lwd=3 )
#lines(density(mydata$af_out), col=alpha("black", 1), lwd=3)
lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80,bw=0.05), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both,bw=0.05), col=alpha("darkorchid2", 1), lwd=3)
#lines(density(avg_perm), col=alpha("darkgoldenrod3", 1), lwd=3, lty=2)

#abline(v=mean(avg_perm), lty=2, col= "black", lwd=3)
#abline(v=mean(snp.sel_75), lty=2, col= "firebrick3", lwd=3)
#abline(v=mean(snp.sel_80), lty=2, col= "royalblue3", lwd=3)
#abline(v=mean(snp.sel_both), lty=2, col= "darkorchid2", lwd=3)

legend("topright", c("permuted", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=1,
    cex=0.7, lwd=2)

mtext(text="B",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)
dev.off()

