# supp figs:

# S2- violin plot
# S3- unfolded af plot
# S5- maf vs af change

#


# this geneartes folded and polarized maf plots. Fig 4 and fig S? respectively.
# also generates the violin plot for supp materials.


library(dplyr)
library(broom)
library(tidyr)
library(scales)
library(gridBase)
library(grid)
library(ggplot2)
library(scales)

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)
cut_off <- 0.05/(9828)

mydata$D1_8_maf <- (sapply(mydata$D1_8_af,function(x)
          ifelse(x > 0.5, (1-x), x)))

selected_7 <- mydata[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off),]
selected_8 <- mydata[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off),]
selected_both <- mydata[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off),]
neutral <- mydata[which(mydata$pH8_selection_pval > cut_off & mydata$pH7_selection_pval > cut_off),]

# pull out allele freqs
selected_7_D7_8 <- selected_7$D7_8_af
selected_7_D7_7 <- selected_7$D7_7_af
selected_7_D1_8 <- selected_7$D1_8_af

selected_8_D7_8 <- selected_8$D7_8_af
selected_8_D7_7 <- selected_8$D7_7_af
selected_8_D1_8 <- selected_8$D1_8_af

selected_both_D7_8 <- selected_both$D7_8_af
selected_both_D7_7 <- selected_both$D7_7_af
selected_both_D1_8 <- selected_both$D1_8_af

selected_all <- unique(rbind(selected_7,selected_8,selected_both))
selected_all_D7_8 <- selected_all$D7_8_af
selected_all_D7_7 <- selected_all$D7_7_af
selected_all_D1_8 <- selected_all$D1_8_af

d7_8 <- (selected_all_D7_8 - selected_all_D1_8)
d7_7 <- (selected_all_D7_7 - selected_all_D1_8)

d7_8_both <- (selected_both_D7_8 - selected_both_D1_8)
d7_7_both <- (selected_both_D7_7 - selected_both_D1_8)

d7_8_s7 <- (selected_7_D7_8 - selected_7_D1_8)
d7_7_s7 <- (selected_7_D7_7 - selected_7_D1_8)

d7_8_s8 <- (selected_8_D7_8 - selected_8_D1_8)
d7_7_s8 <- (selected_8_D7_7 - selected_8_D1_8)


cor.test(d7_7, d7_8,
         method = "pearson",
         conf.level = 0.95)
cor.test(d7_7_s7, d7_8_s7,
         method = "pearson",
         conf.level = 0.95)
cor.test(d7_7_s8, d7_8_s8,
         method = "pearson",
         conf.level = 0.95)

cor.test(d7_7_both, d7_8_both,
         method = "pearson",
         conf.level = 0.95)


# corr genome wide:

d7_8_all <- (neutral$D7_8_af - neutral$D1_8_af)
d7_7_all <- (neutral$D7_7_af - neutral$D1_8_af)

cor.test(d7_8_all, d7_7_all)


######
##
## maf plot
##
######
#sampling from entire dist

#mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)

# need to pull out only selected alleles
snp.sel_75 <- mydata$D1_8_maf[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off)]
snp.sel_80 <- mydata$D1_8_maf[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off)]
snp.sel_both <- mydata$D1_8_maf[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off)]

mean(snp.sel_75)
mean(snp.sel_80)
mean(snp.sel_both)
mean(mydata$D1_8_maf)

sd(snp.sel_75)/sqrt(length(snp.sel_75))
sd(snp.sel_80)/sqrt(length(snp.sel_80))
sd(snp.sel_both)/sqrt(length(snp.sel_both))
sd(mydata$D1_8_maf)/sqrt(length(mydata$D1_8_maf))

nrep <- 1000
nsamp <- length(snp.sel_75)

bs <- matrix(nrow=nsamp, ncol=nrep)
avg_rep <- c()
for(i in 1:nrep){
    bs[,i] <- sample(mydata$D1_8_maf, size=nsamp, replace=FALSE)
    avg_rep <- c(avg_rep, bs[,i])
}

# calc ks tests

p_75_ks <- c()
p_80_ks <- c()
p_both_ks <- c()

for (i in 1:ncol(bs)){
   p_75_ks[i] <- ks.test(snp.sel_75, bs[,i])$p.value
   p_80_ks[i] <- ks.test(snp.sel_80, bs[,i])$p.value
   p_both_ks[i] <- ks.test(snp.sel_both, bs[,i])$p.value
}

ks.test(snp.sel_75, mydata$D1_8_maf)
ks.test(snp.sel_80, mydata$D1_8_maf)
ks.test(snp.sel_both, mydata$D1_8_maf)
ks.test(snp.sel_both, snp.sel_80)
ks.test(snp.sel_both, snp.sel_75)
ks.test(snp.sel_75, snp.sel_80)

length(which(p_75_ks < 0.05))
length(which(p_80_ks < 0.05))
length(which(p_both_ks < 0.05))

# all are 1000

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

densities.qtiles$x <- seq(from=-0.16, to=0.64, by=0.01)


avg_perm <- c()
for (i in 1: ncol(bs)){
    avg_perm <- c(avg_perm, bs[,i])
}



#### plot af change against maf
png("~/urchin_af/figures/Fig_S05_maf_vs_afchange.png", height=85, width=85, units="mm", res=300)

par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
plot(0,type='n', xlim=c(-0.03,.53), ylim=c(-0.03,.39),
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")
box(which="plot")
points(x=selected_7$D1_8_maf, y=d7_7_s7, col=alpha("firebrick3", 0.2),
    bg = alpha("firebrick3", 0.2), cex=0.8, pch=21)
points(x=selected_8$D1_8_maf, y=d7_8_s8, col=alpha("royalblue3", 0.2),
    bg = alpha("royalblue3", 0.2),  cex=0.8, pch=22)
points(x=selected_both$D1_8_maf, y=d7_7_both, col=alpha("darkorchid4", 0.7),
    bg = alpha("darkorchid4", 0.7), cex=0.8, pch=24)
points(x=selected_both$D1_8_maf, y=d7_8_both, col=alpha("darkorchid4", 0.7),
    bg = alpha("darkorchid4", 0.7), cex=0.8, pch=24)

abline(h=0, lty=2, col="grey")
abline(v=0, lty=2, col="grey")
axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(ylab=expression(paste(Delta," allele frequency")), line=1.5, cex.lab=0.9)
title(xlab="Folded minor allele frequency", line=1.5, cex.lab=0.9)

legend("bottomright", c("pH 7.5 selected",
                    "pH 8.0 selected",
                    "Overlapping selected"),
    horiz = FALSE, inset = c(0, 0),
    col = c("firebrick3", "royalblue3","darkorchid4"),
    pt.bg = c("firebrick3", "royalblue3","darkorchid4"),
    pt.cex=1, cex=0.6, pch=c(21,22,24),
    bg="white")

dev.off()




#########
#
# unfolded maf
#
#########


mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)
cut_off <- 0.05/(9828)

# need to pull out only selected alleles
snp.sel_75 <- mydata$D1_8_af[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off)]
snp.sel_80 <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off)]
snp.sel_both <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off)]

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




png("~/urchin_af/figures/Fig_S04_unfolded_af.png", height=85, width=85, units="mm", res=300)
par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)


plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Starting allele frequency", line=1.5, cex.lab=0.9)
title(ylab="Density", line=1.5, cex.lab=0.9)

polygon(x=c(densities.qtiles$x,rev(densities.qtiles$x)),
    y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
    col=alpha("black", alpha=0.2),border=NA)

lines(densities.qtiles$x,densities.qtiles$q50, col="black", lwd=3 )
#lines(density(mydata$af_out), col=alpha("black", 1), lwd=3)
lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80,bw=0.05), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both,bw=0.05), col=alpha("darkorchid2", 1), lwd=3)
abline(v=mean(avg_perm), lty=2, col= "black", lwd=3) # mean: 0.459
abline(v=mean(snp.sel_75)-0.01, lty=2, col= "firebrick3", lwd=3) # mean: 0.248  note, shifting this slightly so visible in plot.
abline(v=mean(snp.sel_80), lty=2, col= "royalblue3", lwd=3) # mean: 0.251
abline(v=mean(snp.sel_both), lty=2, col= "darkorchid2", lwd=3) # mean: 0.32

legend("topright", c("neutral", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=1,
    cex=0.59, lwd=1.9)

dev.off()


# violin plot

d7_8_s8_in <- data.frame(af=d7_8_s8, id=rep("pH 8.0 \nselected", length(d7_8_s8)))
d7_7_s7_in <- data.frame(af=d7_7_s7, id=rep("pH 7.5 \nselected", length(d7_7_s7)))
d7_8_both_in <- data.frame(af=(d7_8_both), id=rep("pH 8.0 \noverlapping", length(d7_8_both)))
d7_7_both_in <- data.frame(af=(d7_7_both), id=rep("pH 7.5 \noverlapping", length(d7_7_both)))
d7_7_neut <- data.frame(af=abs(d7_7_all), id=rep("pH 7.5 \nneutral", length(d7_7_all)))
d7_8_neut <- data.frame(af=abs(d7_8_all), id=rep("pH 8.0 \nneutral", length(d7_8_all)))


new <- rbind(d7_8_neut,d7_8_s8_in,d7_8_both_in, d7_7_neut, d7_7_s7_in, d7_7_both_in)
levels(new$id) <- gsub("overlapping ", "overlapping\n",levels(new$id))

png("~/urchin_af/figures/Fig_S02_violin.png", height=95, width=95, units="mm", res=300)
par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
ggplot(new, aes(factor(id), af, fill=id) ) +
#geom_jitter(height = 0, width = 0.2, alpha = 0.4, size=1)+
geom_violin(alpha = 0.8, draw_quantiles = c(0.5))+
theme_bw() +
theme(legend.position ="none")+
xlab("") + ylab("Change in allele frequency")+
scale_fill_manual(values=c("gray40","royalblue3","darkorchid4",
                    "gray40", "firebrick3", "darkorchid4"))+
#scale_color_manual(values=c("royalblue3", "firebrick3","royalblue3", "firebrick3", "darkorchid4", "darkorchid4"))+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
theme(axis.title.y = element_text(size = rel(0.9)))+
theme(plot.margin = unit(c(5.5,5.5,0.5,5.5), "pt"))

dev.off()

