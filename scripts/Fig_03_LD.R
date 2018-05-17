# Fig_03_LD.R
# Fig S01

library(scales)
library(plyr)
library(dplyr)
library(broom)
library(tidyr)

ld.d1 <- read.table("~/urchin_af/analysis/D1.ld.out", header=FALSE)
ld.d7_7 <- read.table("~/urchin_af/analysis/D7_7.ld.out", header=FALSE)
ld.d7_8 <- read.table("~/urchin_af/analysis/D7_8.ld.out", header=FALSE)

names(ld.d1) <- c("chr", "snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
colnames(ld.d7_7) <-colnames(ld.d1)
colnames(ld.d7_8) <-colnames(ld.d1)

#make id columns for snps
ld.d1$id1 <- paste(ld.d1$chr, ld.d1$snp1, sep=":")
ld.d1$id2 <- paste(ld.d1$chr, ld.d1$snp2, sep=":")

ld.d7_7$id1 <- paste(ld.d7_7$chr, ld.d7_7$snp1, sep=":")
ld.d7_7$id2 <- paste(ld.d7_7$chr, ld.d7_7$snp2, sep=":")

ld.d7_8$id1 <- paste(ld.d7_8$chr, ld.d7_8$snp1, sep=":")
ld.d7_8$id2 <- paste(ld.d7_8$chr, ld.d7_8$snp2, sep=":")

#read in cmh results
cmh <-read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

cmh$sel_pH7 <- FALSE
cmh$sel_pH8<- FALSE

# thin to those that are showing response in treatment
cmh$sel_pH7[which(cmh$pH7_selection_qval < 0.001)] <- TRUE

cmh$sel_pH8[which(cmh$pH8_selection_qval < 0.001)] <- TRUE

cmh$SNP <- paste(cmh$CHROM, cmh$POS, sep=":")

selected <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, cmh$sel_pH7, cmh$sel_pH8 )
colnames(selected) <- c("chr", "pos", "id", "sel_pH7", "sel_pH8")

dat.d1 <- merge(ld.d1, selected, by.x="id1", by.y="id", all.x=TRUE)
dat.d1 <- merge(dat.d1, selected, by.x="id2", by.y="id", all.x=TRUE)

dat.d7_7 <- merge(ld.d7_7, selected, by.x="id1", by.y="id", all.x=TRUE)
dat.d7_7 <- merge(dat.d7_7, selected, by.x="id2", by.y="id", all.x=TRUE)

dat.d7_8 <- merge(ld.d7_8, selected, by.x="id1", by.y="id", all.x=TRUE)
dat.d7_8 <- merge(dat.d7_8, selected, by.x="id2", by.y="id", all.x=TRUE)

##############
###
### plot decay in ld
###
##############

# calc distance between snp1 and snp2
dat.d1$distance <- abs(dat.d1$snp1-dat.d1$snp2)
dat.d7_7$distance <- abs(dat.d7_7$snp1-dat.d7_7$snp2)
dat.d7_8$distance <- abs(dat.d7_8$snp1-dat.d7_8$snp2)

dat.d1.1 <- dat.d1
dat.d7_7.1 <- dat.d7_7
dat.d7_8.1 <- dat.d7_8

dat.d7_7 <- dat.d7_7[which(dat.d7_7$distance <=200),]
dat.d7_8 <- dat.d7_8[which(dat.d7_8$distance <=200),]
dat.d1 <- dat.d1[which(dat.d1$distance <=200),]

dat.d7_7_sel <- dat.d7_7[which(dat.d7_7$sel_pH7.y == TRUE | dat.d7_7$sel_pH7.x == TRUE),]
dat.d7_7_neut <- dat.d7_7[which(dat.d7_7$sel_pH7.y == FALSE & dat.d7_7$sel_pH7.x == FALSE),]

dat.d7_8_sel <- dat.d7_8[which(dat.d7_8$sel_pH8.y == TRUE | dat.d7_8$sel_pH8.x == TRUE),]
dat.d7_8_neut <- dat.d7_8[which(dat.d7_8$sel_pH8.y == FALSE & dat.d7_8$sel_pH8.x == FALSE),]

mod.d1 <- nls(mle_est ~ exp(a + b * distance), data = dat.d1, start = list(a = 0, b = 0))
mod.d7_7 <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_7, start = list(a = 0, b = 0))
mod.d7_8 <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_8, start = list(a = 0, b = 0))

#plot low pH selected vs neutral
mod.d7_7_sel <-  nls(mle_est ~ exp(a + b * distance), data = dat.d7_7_sel,  start = list(a = 0, b = 0))
mod.d7_7_neut <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_7_neut, start = list(a = 0, b = 0))
mod.d7_8_sel <-  nls(mle_est ~ exp(a + b * distance), data = dat.d7_8_sel,  start = list(a = 0, b = 0))
mod.d7_8_neut <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_8_neut, start = list(a = 0, b = 0))

## sample random subgroups of data. compare to selected loci

perm_length <- length(which(cmh$sel_pH7 == TRUE))

ks.val.sel <- c()
ks.mod.sel <- c()
ks.val.ctr <- c()
ks.mod.ctr <- c()

perm.out <- list()

for (i in 1:500){
    rand.sel <- rep("FALSE", nrow(cmh))
    rand.sel[sample(1:nrow(cmh), perm_length)]<- TRUE
    selected.perm <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, rand.sel )
    colnames(selected.perm) <- c("chr", "pos", "id", "selected")
    dat.perm <- merge(ld.d7_7, selected.perm, by.x="id1", by.y="id", all.x=TRUE)
    dat.perm <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
    dat.perm$distance <- abs(dat.perm$snp1-dat.perm$snp2)
    dat.perm <- dat.perm[which(dat.perm$distance <=200),]
    dat.perm <- dat.perm[which(dat.perm$selected.y == TRUE | dat.perm$selected.x == TRUE),]
    mod.perm <- nls(mle_est ~ exp(a + b * distance), data = dat.perm, start = list(a = 0, b = 0))
    #lines(sort(dat.perm$distance, decreasing=FALSE),
    #    sort(fitted(mod.perm), decreasing=TRUE),
    #    lwd=1, col=rgb(0,0,0,alpha=0.1))
    if (i%%10 == 0){print(i)} # printing progress
    ks.mod.sel[i] <- ks.test(fitted(mod.d7_7_sel), fitted(mod.perm))$p.value
    ks.val.sel[i] <- ks.test(dat.d7_7_sel$mle_est, dat.perm$mle_est)$p.value
    ks.mod.ctr[i] <- ks.test(fitted(mod.d7_8_sel), fitted(mod.perm))$p.value
    ks.val.ctr[i] <- ks.test(dat.d7_8_sel$mle_est, dat.perm$mle_est)$p.value
    perm.out[[i]] <- cbind(sort(dat.perm$distance, decreasing=FALSE), sort(fitted(mod.perm), decreasing=TRUE))

}


#########
##
## plot actual fig 2
##
#########
## figure out CI

# permutations saved in: perm.out

# for each permutation, take the mean of the bin, then combine all perms.

out <- list()
for (i in 1: length(perm.out)){

    tmp.df <- data.frame(x=perm.out[[i]][,1],
                            y= perm.out[[i]][,2],
                        bin = cut(
                            perm.out[[i]][,1],
                            breaks=seq(from=0, to=200, by=1),
                            labels=seq(from=1, to=200, by=1)))
    out[[i]] <- data.frame(bin=seq(from=1, to=200, by=1), ld=tapply(tmp.df$y,tmp.df$bin , mean))
}

out.new <- do.call(rbind, out)

densities.qtiles <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(ld, 0.025, na.rm=TRUE),
            q50 = quantile(ld, 0.5, na.rm=TRUE),
            q95 = quantile(ld, 0.975, na.rm=TRUE))


tiff("~/urchin_af/figures/Fig_03_ld.tiff", height=85, width=85, units="mm", res=300)
par(mfrow = c(1, 1), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(0,type='n', xlim=c(0,200), ylim=c(0,0.4),
    main="",
    ylab="",
    xlab="",
    cex.lab=0.9, cex.axis=0.7,
    xaxt="n",yaxt="n")

polygon(x=c(densities.qtiles$bin,rev(densities.qtiles$bin)),
    y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
    col=alpha("black", alpha=0.3),border=NA)

lines(densities.qtiles$bin, densities.qtiles$q50, lwd=3, lty=2, col='black')
lines(sort(dat.d7_8_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut), decreasing=TRUE),
    lwd=3.2, lty=2, col='royalblue3')
lines(sort(dat.d7_7_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut), decreasing=TRUE),
    lwd=3, lty=2, col='firebrick3')

# selected sites
lines(sort(dat.d7_8_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel), decreasing=TRUE),
    lwd=3, col='royalblue3')
lines(sort(dat.d7_7_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel), decreasing=TRUE),
    lwd=3, col='firebrick3')

axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)
title(xlab="Distance between SNPs in base pairs", line=2, cex.lab=0.9)
title(ylab="Estimated Linkage Disequilibrium", line=2, cex.lab=0.9)

legend("topright",
    legend=c(expression('T'[0]),
            "pH 7.5: selected",
            "pH 7.5: neutral",
            "pH 8.0: selected",
            "pH 8.0: neutral"),
       col=c("black", "firebrick3", "firebrick3", "royalblue3", "royalblue3"),
       lty=c(2,1,2,1,2), lwd=2.2, cex=0.8)

dev.off()




png("~/urchin_af/figures/Fig_S02_ldperm.png", height=100, width=200, units="mm", res=300)
#par(mfrow = c(1, 2), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)


plot(0,type='n', xlim=c(0,200), ylim=c(0,0.4),
    main="",
    ylab="",
    xlab="",
    cex.lab=0.9, cex.axis=0.7,
    xaxt="n",yaxt="n")


for (i in 1:length(perm.out)){
    lines(perm.out[[i]][,1],
          perm.out[[i]][,2],
        lwd=1, col=rgb(0,0,0,alpha=0.1))
   if (i%%25 == 0){print(i)} # printing progress
}


lines(sort(dat.d7_7_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut), decreasing=TRUE),
    lwd=3.1, lty=2, col='firebrick3')
lines(sort(dat.d7_8_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut), decreasing=TRUE),
    lwd=2.9, lty=2, col='royalblue3')
lines(sort(dat.d7_8_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel), decreasing=TRUE),
    lwd=3, col='royalblue3')
lines(sort(dat.d7_7_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel), decreasing=TRUE),
    lwd=3, col='firebrick3')


legend("topright", c(expression('T'[0] ~ "Permuted"),
            "pH 7.5: selected",
            "pH 7.5: neutral",
            "pH 8.0: selected",
            "pH 8.0: neutral"),
       col=c("black", "firebrick3", "firebrick3", "royalblue3", "royalblue3"),
       lty=c(1,2,1,2,1), lwd=2.2, cex=0.8)


dev.off()




##############################################################################
#### plot by bin avg:
##############################################################################

#revert to all distances
dat.d7_7 <- dat.d7_7.1
dat.d7_8 <- dat.d7_8.1
dat.d1 <- dat.d1.1

dat.d7_7_sel <- dat.d7_7[which(dat.d7_7$sel_pH7.y == TRUE | dat.d7_7$sel_pH7.x == TRUE),]
dat.d7_7_neut <- dat.d7_7[which(dat.d7_7$sel_pH7.y == FALSE & dat.d7_7$sel_pH7.x == FALSE),]

dat.d7_8_sel <- dat.d7_8[which(dat.d7_8$sel_pH8.y == TRUE | dat.d7_8$sel_pH8.x == TRUE),]
dat.d7_8_neut <- dat.d7_8[which(dat.d7_8$sel_pH8.y == FALSE & dat.d7_8$sel_pH8.x == FALSE),]

dat.d7_7_sel$bin <- cut(dat.d7_7_sel$distance, breaks=seq(from=0,to=500, by=10))
dat.d7_7_sel$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_7_sel$bin),","), '[', 1))
dat.d7_7_neut$bin <- cut(dat.d7_7_neut$distance, breaks=seq(from=0,to=500, by=10))
dat.d7_7_neut$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_7_neut$bin),","), '[', 1))

dat.d7_8_sel$bin <- cut(dat.d7_8_sel$distance, breaks=seq(from=0,to=500, by=10))
dat.d7_8_sel$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_8_sel$bin),","), '[', 1))
dat.d7_8_neut$bin <- cut(dat.d7_8_neut$distance, breaks=seq(from=0,to=500, by=10))
dat.d7_8_neut$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_8_neut$bin),","), '[', 1))

dat.d1$bin <- cut(dat.d1$distance, breaks=seq(from=0,to=500, by=10))
dat.d1$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d1$bin),","), '[', 1))

a <- aggregate(mle_est ~ bin1 , data=dat.d7_7_sel, mean)
a.ct <- aggregate(mle_est ~ bin1 , data=dat.d7_7_sel, length)
a <- data.frame(bin1=a$bin1, mle_est=a$mle_est, count= a.ct$mle_est)

b <- aggregate(mle_est ~ bin1 , data=dat.d7_7_neut, mean)
b.ct <- aggregate(mle_est ~ bin1 , data=dat.d7_7_neut, length)
b <- data.frame(bin1=b$bin1, mle_est=b$mle_est, count= b.ct$mle_est)

c <- aggregate(mle_est ~ bin1 , data=dat.d7_8_sel, mean)
c.ct <- aggregate(mle_est ~ bin1 , data=dat.d7_8_sel, length)
c <- data.frame(bin1=c$bin1, mle_est=c$mle_est, count= c.ct$mle_est)

d <- aggregate(mle_est ~ bin1 , data=dat.d7_8_neut, mean)
d.ct <- aggregate(mle_est ~ bin1 , data=dat.d7_8_neut, length)
d <- data.frame(bin1=d$bin1, mle_est=d$mle_est, count= d.ct$mle_est)

e <- aggregate(mle_est ~ bin1 , data=dat.d1, mean)
e.ct <- aggregate(mle_est ~ bin1 , data=dat.d1, length)
e <- data.frame(bin1=e$bin1, mle_est=e$mle_est, count= e.ct$mle_est)

# with transparency scaled by number of variants

range01 <- function(x){
    if(max(x)>= 2.69897){ # log10(500) for cutoff. this could be changed
        a <- x
        a[which(a > 2.69897)] <- 2.69897
        out <- (a-min(a))/((max(a)-min(a)))
        return(out)}
    else{(x-min(x))/((max(x)-min(x)))}
    }


a$count.scale <- range01(log10(a$count))
b$count.scale <- range01(log10(b$count))
c$count.scale <- range01(log10(c$count))
d$count.scale <- range01(log10(d$count))
e$count.scale <- range01(log10(e$count))

png("~/urchin_af/figures/Fig_S01_ldbinned.png", res=300, height=4, width=7, units="in")
par(mfrow = c(1, 1),mar=c(4.5, 4.5, 0.5, 0.5))

plot(0,type='n', xlim=c(0,500), ylim=c(0,0.5),
    xlab="Distance between SNPs", ylab="Binned average estimated LD",
    main="", cex.axis=0.8, cex.lab=0.9 )

for(i in 1:nrow(e)){
    points(x=as.numeric(as.character(e$bin1[i])), y=as.numeric(as.character(e$mle_est[i])),
        col="black", bg=rgb(0,0,0,alpha=e$count.scale[i]),
        pch=22, cex=1.1, lwd=1.5)
}

for(i in 1:nrow(b)){
    points(x=as.numeric(as.character(b$bin1[i])), y=as.numeric(as.character(b$mle_est[i])),
        col=rgb(0.75,0,0,alpha=1), bg=rgb(0.75,0,0,alpha=b$count.scale[i]),
        pch=21, cex=1.1, lwd=1.5)
}

for(i in 1:nrow(d)){
    points(x=as.numeric(as.character(d$bin1[i])), y=as.numeric(as.character(d$mle_est[i])),
        col=rgb(0,0,0.75,alpha=1), bg=rgb(0,0,0.75,alpha=d$count.scale[i]),
        pch=22, cex=1.1, lwd=1.5)
}

for(i in 1:nrow(c)){
    points(x=as.numeric(as.character(c$bin1[i])), y=as.numeric(as.character(c$mle_est[i])),
        col=rgb(0,0,0.75,alpha=1), bg=rgb(0,0,0.75,alpha=c$count.scale[i]),
        pch=24, cex=1.1, lwd=1.5)
}

for(i in 1:nrow(a)){
    points(x=as.numeric(as.character(a$bin1[i])), y=as.numeric(as.character(a$mle_est[i])),
        col=rgb(0.75,0,0,alpha=1), bg=rgb(0.75,0,0,alpha=a$count.scale[i]),
        pch=25, cex=1.1, lwd=1.5)
}

legend(75,0.5,legend=c("pH 7.5: selected","pH 7.5: neutral",
    "pH 8.0: selected", "pH 8.0: neutral", expression('T'[0])),
       col=c("red", "red", "blue","blue", "black"), pt.bg=c("red", "red", "blue","blue", "black"),
        pch=c(25,21,24,22,22), cex=0.9)


dev.off()
