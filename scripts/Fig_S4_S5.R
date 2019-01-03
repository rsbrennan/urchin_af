### Fig S5, S6

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
perm_dat <- read.table("~/urchin_af/analysis/permutation_af.txt", stringsAsFactors=FALSE, header=TRUE)
cut_off <- (0.05/9828)
# need to pull out only selected alleles
snp.sel_75 <- mydata$D1_8_af[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off)]
snp.sel_80 <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off)]
snp.sel_both <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off)]

######
##
## Fig S4
##
######

# read in data from permutation in 06_balancing_sel_permutation.R

## plot CI

dat.bs <- data.frame(CHROM = mydata$CHROM, POS = mydata$POS, af_out = mydata$D1_8_af)

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

png("~/urchin_af/figures/Fig_S4_cmhperm.png", height=100, width=100, units="mm", res=300)

par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
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


# not longer using, individual lines from perm:


#plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
#    main="",
#    ylab="",
#    xlab="",
#    cex.lab=1.1, cex.axis=1,
#    xaxt="n",yaxt="n")
#
#axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
#axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
#title(ylab="Density", line=1.5, cex.lab=1)
#title(xlab="Starting allele frequency", line=1.5, cex.lab=1)
#
#avg_perm <- c()
#for (i in 1: ncol(perm_dat)){
#    lines(density(unlist(perm_dat[,i]),na.rm=TRUE, bw=0.05), col=alpha("black", 0.08), lwd=1)
#    avg_perm <- c(avg_perm, unlist(perm_dat[,i]))
#}
#
#lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
#lines(density(snp.sel_80, bw=0.05), col=alpha("royalblue3", 1), lwd=3)
#lines(density(snp.sel_both, bw=0.05), col=alpha("darkorchid2", 1), lwd=3)
##lines(density(avg_perm), col=alpha("darkgoldenrod3", 1), lwd=3, lty=2)
#
##abline(v=mean(avg_perm), lty=2, col= "black", lwd=3)
##abline(v=mean(snp.sel_75), lty=2, col= "firebrick3", lwd=3)
##abline(v=mean(snp.sel_80), lty=2, col= "royalblue3", lwd=3)
##abline(v=mean(snp.sel_both), lty=2, col= "darkorchid2", lwd=3)
#
#legend("topright", c("permuted", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
#    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=1,
#    cex=0.8, lwd=2)
#
#
#mtext(text="A",
#        side=3, line=0,
#             cex=1.5,
#            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)
#

##################################################
##### Fig S5
####################################################### 

# all lines of samps from genome wide maf
# hist of MAF

all_bin <- cut(mydata$D1_8_af, breaks=seq(0,1, 0.05))
snp.sel_80_bin <- cut(snp.sel_80, breaks=seq(0,1, 0.05))
snp.sel_75_bin <- cut(snp.sel_75, breaks=seq(0,1, 0.05))
snp.sel_both_bin <- cut(snp.sel_both, breaks=seq(0,1, 0.05))

bin <- c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25",
    "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5",
    "0.5-0.55", "0.55-0.6", "0.6-0.65", "0.65-0.7", "0.7-0.75",
    "0.75-0.8", "0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-1")

# make barplot of neutral vs selected
all_freq <- table(all_bin)/sum(table(all_bin))
snp.sel_80_freq <- table(snp.sel_80_bin)/sum(table(snp.sel_80_bin))
snp.sel_75_freq <- table(snp.sel_75_bin)/sum(table(snp.sel_75_bin))
snp.sel_both_freq <- table(snp.sel_both_bin)/sum(table(snp.sel_both_bin))

# merge into tables for plotting
all_freq.dat <- data.frame(class=rep("all", length(bin)),bin=bin,
            allele_frequency=as.vector(all_freq))
snp.sel_80.dat <- data.frame(class=rep("pH 8.0", length(bin)),bin=bin,
            allele_frequency=as.vector(snp.sel_80_freq))
snp.sel_75.dat <- data.frame(class=rep("pH 7.5", length(bin)),bin=bin,
            allele_frequency=as.vector(snp.sel_75_freq))
snp.sel_both.dat <- data.frame(class=rep("Overlapping selected", length(bin)),bin=bin,
            allele_frequency=as.vector(snp.sel_both_freq))

fq.dat <- rbind(all_freq.dat,snp.sel_80.dat, snp.sel_75.dat, snp.sel_both.dat)

png("~/urchin_af/figures/Fig_S5_mafcomp.png", height=100, width=200, units="mm", res=300)

plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
par(fig = gridFIG(), new = TRUE, mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)


#Draw ggplot
p <- ggplot(data=fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('grey38','royalblue3','firebrick3', 'darkorchid2'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.y = element_text(size = rel(0.9)))+
  ylim(0, 0.45) +
  guides(size = FALSE) +
  ylab("relative frequency") +
  theme(legend.title=element_blank())+
  theme(legend.position="none")

print(p, newpage = FALSE)
popViewport()

####################
#### B
####################
# all lines of samps from genome wide maf

#Draw bsae plot
pushViewport(viewport(layout.pos.col = 2))
par(fig = gridFIG(), new = TRUE,mar=c(4, 3, 0.5, 1), mgp=c(3, 1, 0), las=0)

#par(mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Starting allele frequency", line=1.5, cex.lab=1)
title(ylab="Density", line=1.5, cex.lab=1)

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", stringsAsFactors=FALSE, header=TRUE)

nrep <- 1000
nsamp <- length(snp.sel_75)

bs <- matrix(nrow=nsamp, ncol=nrep)
avg_rep <- c()
for(i in 1:nrep){
    bs[,i] <- sample(mydata$D1_8_af, size=nsamp, replace=FALSE)
    avg_rep <- c(avg_rep, bs[,i])
}

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


for (i in 1: ncol(bs)){
    lines(density(bs[,i], bw=0.05), col=alpha("black", 0.08), lwd=1)
}
#lines(density(mydata$af_out), col=alpha("black", 1), lwd=3)
lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80,bw=0.05), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both,bw=0.05), col=alpha("darkorchid2", 1), lwd=3)

legend("topright", c("All variants", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), 
    lty=1,
    cex=0.8, lwd=2)

popViewport()

dev.off()
