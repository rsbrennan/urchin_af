library(scales)

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

cut_off <- quantile(cmh$control_selection_pval, 0.01, na.rm=TRUE)
cut_d7 <- quantile(cmh$d7_selection_pval, 0.01, na.rm=TRUE)

cut_off
cmh$selected <- FALSE
cmh$d7_sel<- FALSE
cmh$lab_sel <- FALSE
# thin to those that are showing response in treatment
cmh$selected[which(cmh$pH_selection_pval < cut_off)] <- TRUE
#cmh$lab_sel[which(cmh$control_selection_pval < cut_off)] <- TRUE
#cmh$d7_sel[which(cmh$d7_selection_pval < cut_d7)] <- TRUE
cmh$d7_sel[which(cmh$d7_selection_pval < cut_off)] <- TRUE

cmh$lab_sel[sample(1:nrow(cmh), length(which(cmh$pH_selection_pval < cut_off)))] <- TRUE


cmh$SNP <- paste(cmh$CHROM, cmh$POS, sep=":")

selected <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, cmh$selected, cmh$lab_sel, cmh$d7_sel )
colnames(selected) <- c("chr", "pos", "id", "selected", "lab_sel", "d7_sel")

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

dat.d7_7 <- dat.d7_7[which(dat.d7_7$distance <=200),]
dat.d7_8 <- dat.d7_8[which(dat.d7_8$distance <=200),]
dat.d1 <- dat.d1[which(dat.d1$distance <=200),]

dat.d7_7_sel <- dat.d7_7[which(dat.d7_7$selected.y == TRUE | dat.d7_7$selected.x == TRUE),]
dat.d7_7_neut <- dat.d7_7[which(dat.d7_7$selected.y == FALSE & dat.d7_7$selected.x == FALSE),]

dat.d7_d7_sel <- dat.d7_7[which(dat.d7_7$d7_sel.y == TRUE | dat.d7_7$d7_sel.x == TRUE),]
dat.d7_d7_neut <- dat.d7_7[which(dat.d7_7$d7_sel.y == FALSE & dat.d7_7$d7_sel.x == FALSE),]

dat.d7_8_sel <- dat.d7_8[which(dat.d7_8$lab_sel.y == TRUE | dat.d7_8$lab_sel.x == TRUE),]
dat.d7_8_neut <- dat.d7_8[which(dat.d7_8$lab_sel.y == FALSE & dat.d7_8$lab_sel.x == FALSE),]


#plot(x=dat.d1$distance, y=dat.d1$mle_est)

mod.d1 <- nls(mle_est ~ exp(a + b * distance), data = dat.d1, start = list(a = 0, b = 0))
mod.d7_7 <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_7, start = list(a = 0, b = 0))
mod.d7_8 <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_8, start = list(a = 0, b = 0))
mod.d7_d7 <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_7, start = list(a = 0, b = 0))

#plot low pH selected vs neutral
mod.d7_7_sel <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_7_sel, start = list(a = 0, b = 0))
mod.d7_7_neut <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_7_neut, start = list(a = 0, b = 0))
mod.d7_8_sel <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_8_sel, start = list(a = 0, b = 0))
mod.d7_8_neut <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_8_neut, start = list(a = 0, b = 0))

mod.d7_d7_sel <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_d7_sel, start = list(a = 0, b = 0))
mod.d7_d7_neut <- nls(mle_est ~ exp(a + b * distance), data = dat.d7_d7_neut, start = list(a = 0, b = 0))

png("~/urchin_af/figures/ld_decay_4.png", res=300, height=7, width=7, units="in")
plot(0,type='n', xlim=c(0,200), ylim=c(0,0.4),
    xlab="Distance between SNPs", ylab="Estimated LD",
    main="Decay in LD with distance")

# add day 1
lines(sort(dat.d1$distance, decreasing=FALSE), sort(fitted(mod.d1), decreasing=TRUE), lwd=3, lty=2, col='black')

lines(sort(dat.d7_8_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut), decreasing=TRUE), lwd=3, lty=2, col='blue')
lines(sort(dat.d7_7_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut), decreasing=TRUE), lwd=3, lty=2, col='red')
lines(sort(dat.d7_d7_neut$distance, decreasing=FALSE), sort(fitted(mod.d7_d7_neut), decreasing=TRUE), lwd=3, lty=2, col='orange')

lines(sort(dat.d7_8_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel), decreasing=TRUE), lwd=3, col='blue')

lines(sort(dat.d7_7_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel), decreasing=TRUE), lwd=3, col='red')
lines(sort(dat.d7_d7_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_d7_sel), decreasing=TRUE), lwd=3, col='orange')


legend(80,0.4,
    legend=c("pH day 7: selected","pH day 7: neutral", 
            "Control day 7: selected", "Control day 7: neutral", "Control day 1: neutral",
            "D7 selected", "D7 neutral"),
       col=c("red", "red", "blue", "blue", "black", "orange", "orange"), lty=c(1,2,1,2,2, 1, 2), lwd=3, cex=1.1)

dev.off()

## sample random subgroups of data. compare to selected loci

png("~/urchin_af/figures/ld_1.png", res=300, height=7, width=10, units="in")

par(mfrow = c(1, 2))

plot(0,type='n', xlim=c(0,200), ylim=c(0,0.4),
    xlab="Distance between SNPs", ylab="Estimated LD",
    main="Decay in LD with distance")
a <- c()
b <- c()
ks.val <- c()
ks.mod <- c()
for (i in 1:500){
    rand.sel <- rep("FALSE", nrow(cmh))
    rand.sel[sample(1:nrow(cmh), length(which(cmh$pH_selection_pval < cut_off)))]<- TRUE
    selected.perm <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, rand.sel )
    colnames(selected.perm) <- c("chr", "pos", "id", "selected")
    dat.perm <- merge(ld.d7_7, selected.perm, by.x="id1", by.y="id", all.x=TRUE)
    dat.perm <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
    dat.perm$distance <- abs(dat.perm$snp1-dat.perm$snp2)
    dat.perm <- dat.perm[which(dat.perm$distance <=200),]
    dat.perm <- dat.perm[which(dat.perm$selected.y == TRUE | dat.perm$selected.x == TRUE),]
    mod.perm <- nls(mle_est ~ exp(a + b * distance), data = dat.perm, start = list(a = 0, b = 0))
    lines(sort(dat.perm$distance, decreasing=FALSE), 
        sort(fitted(mod.perm), decreasing=TRUE), 
        lwd=1, col=rgb(0,0,0,alpha=0.2))
   if (i%%10 == 0){print(i)} # printing progress
    a[i] <- coef(mod.perm)["a"]
    b[i] <- coef(mod.perm)["b"]
    ks.mod[i] <- ks.test(fitted(mod.d7_7_sel), fitted(mod.perm))$p.value
    ks.val[i] <- ks.test(dat.d7_7_sel$mle_est, dat.perm$mle_est)$p.value

}
    
lines(sort(dat.d7_7_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel), decreasing=TRUE), lwd=3, col='red')

plot(0,type='n', xlim=c(0,200), ylim=c(0,0.4),
    xlab="Distance between SNPs", ylab="Estimated LD",
    main="Decay in LD with distance")

c <- c()
d <- c()
ks.val.a <- c()
ks.mod.a <- c()
for (i in 1:500){
    rand.sel <- rep("FALSE", nrow(cmh))
    rand.sel[sample(1:nrow(cmh), length(which(cmh$d7_selection_pval < cut_d7)))]<- TRUE
    selected.perm <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, rand.sel )
    colnames(selected.perm) <- c("chr", "pos", "id", "selected")
    dat.perm <- merge(ld.d7_7, selected.perm, by.x="id1", by.y="id", all.x=TRUE)
    dat.perm <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
    dat.perm$distance <- abs(dat.perm$snp1-dat.perm$snp2)
    dat.perm <- dat.perm[which(dat.perm$distance <=200),]
    dat.perm <- dat.perm[which(dat.perm$selected.y == TRUE | dat.perm$selected.x == TRUE),]
    mod.perm <- nls(mle_est ~ exp(a + b * distance), data = dat.perm, start = list(a = 0, b = 0))
    lines(sort(dat.perm$distance, decreasing=FALSE), 
        sort(fitted(mod.perm), decreasing=TRUE), 
        lwd=1, col=rgb(0,0,0,alpha=0.2))
    if (i%%10 == 0){print(i)} # printing progress
    c[i] <- coef(mod.perm)["a"]
    d[i] <- coef(mod.perm)["b"]
    ks.mod.a[i] <- ks.test(fitted(mod.d7_d7_sel), fitted(mod.perm))$p.value
    ks.val.a[i] <- ks.test(dat.d7_d7_sel$mle_est, dat.perm$mle_est)$p.value
}

lines(sort(dat.d7_d7_sel$distance, decreasing=FALSE), sort(fitted(mod.d7_d7_sel), decreasing=TRUE), lwd=3, col='orange')

dev.off()   

# plot ks test results
png("~/urchin_af/figures/ld_ks.png", res=300, height=7, width=10, units="in")

par(mfrow = c(2, 2))

hist(ks.val , breaks=40, col=rgb(0,0,0,alpha=0.2), 
    main="Day1 vs Day 7 ks tests:\n raw LD vales")

hist(ks.mod , breaks=40, col=rgb(0,0,0,alpha=0.2), 
    main="Day1 vs Day 7 ks tests:\n vals predicted from LD model")

hist(ks.val.a , breaks=40, col=rgb(0,0,0,alpha=0.2), 
    main="Day 7 comparison ks tests:\n raw LD vales")

hist(ks.mod.a , breaks=40, col=rgb(0,0,0,alpha=0.2), 
    main="Day 7 comparison ks tests:\n vals predicted from LD model")

dev.off()


png("~/urchin_af/figures/ld_3.png", res=300, height=7, width=10, units="in")
par(mfrow = c(2, 2))

p <- coef(mod.d7_7_sel)

hist(a, breaks=20, col='grey', main="D1 vs. D7: hist of a")
abline(v=p["a"], col="red", lwd=2)

hist(b, breaks=20, col='grey', main="D1 vs. D7: hist of b")
abline(v=p["b"], col="red", lwd=2)

p <- coef(mod.d7_d7_sel)
hist(c, breaks=20, col='grey', main="D7 vs. D7: hist  of a")
abline(v=p["a"], col="orange", lwd=2)

hist(d, breaks=20, col='grey', main="D7 vs. D7: hist  of b")
abline(v=p["b"], col="orange", lwd=2)
dev.off()

########################################################################################
############################################
##
## generate/merge outliers
## 
##
############################################
########################################################################################

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")
cmh$SNP <- paste(cmh$CHROM, cmh$POS, sep=":")

cut_off <- quantile(cmh$control_selection_pval, 0.01, na.rm=TRUE)
cut_d7 <- quantile(cmh$d7_selection_pval, 0.01, na.rm=TRUE)

p1 <- mydata$SNP[which(mydata$pH_selection_pval < cut_off)]
#p2 <- mydata$SNP[which(mydata$d7_selection_pval < quantile(mydata$d7_selection_pval, 0.01, na.rm=TRUE))]
p2 <- mydata$SNP[which(mydata$d7_selection_pval < cut_off)]

inter.sel <- intersect(p1, p2)
overlap.sel <- rep("FALSE", nrow(mydata))
overlap.sel[(mydata$SNP %in% inter.sel)]<- TRUE
overlap.sel <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, overlap.sel )
colnames(overlap.sel) <- c("chr", "pos", "id", "selected")
overlap.new <- merge(ld.d7_7, overlap.sel, by.x="id1", by.y="id", all.x=TRUE)
overlap.new <- merge(overlap.new, overlap.sel, by.x="id2", by.y="id", all.x=TRUE)
overlap.new$distance <- abs(overlap.new$snp1-overlap.new$snp2)
overlap.new <- overlap.new[which(overlap.new$distance <=200),]
overlap.new <- overlap.new[which(overlap.new$selected.y == TRUE | overlap.new$selected.x == TRUE),]
mod.overlap <- nls(mle_est ~ exp(a + b * distance), data = overlap.new, start = list(a = 0, b = 0))

png("~/urchin_af/figures/ld_2.png", res=300, height=7, width=10, units="in")
par(mfrow = c(1, 1))

plot(0,type='n', xlim=c(0,200), ylim=c(0,0.6),
    xlab="Distance between SNPs", ylab="Estimated LD",
    main="overlapping variants ld decay")
a <- c()
b <- c()
for (i in 1:300){
    rand.sel <- rep("FALSE", nrow(cmh))
    rand.sel[sample(1:nrow(cmh), length(inter.sel))]<- TRUE
    selected.perm <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, rand.sel )
    colnames(selected.perm) <- c("chr", "pos", "id", "selected")
    dat.perm <- merge(ld.d7_7, selected.perm, by.x="id1", by.y="id", all.x=TRUE)
    dat.perm <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
    dat.perm$distance <- abs(dat.perm$snp1-dat.perm$snp2)
    dat.perm <- dat.perm[which(dat.perm$distance <=200),]
    dat.perm <- dat.perm[which(dat.perm$selected.y == TRUE | dat.perm$selected.x == TRUE),]
    mod.perm <- nls(mle_est ~ exp(a + b * distance), data = dat.perm, start = list(a = 0, b = 0), control = list(maxiter = 500))
    lines(sort(dat.perm$distance, decreasing=FALSE), 
        sort(fitted(mod.perm), decreasing=TRUE), 
        lwd=1, col=rgb(0,0,0,alpha=0.2))
    a[i] <- coef(mod.perm)["a"]
    b[i] <- coef(mod.perm)["b"]
   if (i%%10 == 0){print(i)} # printing progress
}

lines(sort(overlap.new$distance, decreasing=FALSE), sort(fitted(mod.overlap), decreasing=TRUE), lwd=3.5, col='green')
dev.off()

png("~/urchin_af/figures/ld_4.png", res=300, height=7, width=10, units="in")
p <- coef(mod.overlap)
par(mfrow = c(1, 2))

hist(a, breaks=50, col='grey', main="overlap: hist of a", xlim=c(-3.5,1))
abline(v=p["a"], col="green", lwd=2)

hist(b, breaks=60, col='grey', main="overlap: hist of b", xlim=c(-0.1,0.05))
abline(v=p["b"], col="green", lwd=2)

dev.off()

# Save coefficients in a vector

#save model

fx <- function(a,b, distance)exp(a + b * distance)
# Plot and add a curve. Note that 'p' is a named vector.

lines(sort(overlap.new$distance, decreasing=FALSE), 
        fx(a=-1, b=-0.01, distance=sort(overlap.new$distance, decreasing=FALSE)),
        col='red', lwd=2)

# a is the intercept. 
# b is the depth of the curve. where it asympotyes


#############
#### plot by bin avg:
#############

library(plyr)

dat.d7_7_sel$bin <- cut(dat.d7_7_sel$distance, breaks=seq(from=0,to=300, by=10))
dat.d7_7_sel$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_7_sel$bin),","), '[', 1))
dat.d7_7_neut$bin <- cut(dat.d7_7_neut$distance, breaks=seq(from=0,to=300, by=10))
dat.d7_7_neut$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_7_neut$bin),","), '[', 1))


dat.d7_d7_sel$bin <- cut(dat.d7_d7_sel$distance, breaks=seq(from=0,to=300, by=10))
dat.d7_d7_sel$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_d7_sel$bin),","), '[', 1))

dat.d7_d7_neut$bin <- cut(dat.d7_d7_neut$distance, breaks=seq(from=0,to=300, by=10))
dat.d7_d7_neut$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_d7_neut$bin),","), '[', 1))

dat.d7_8_sel$bin <- cut(dat.d7_8_sel$distance, breaks=seq(from=0,to=300, by=10))
dat.d7_8_sel$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_8_sel$bin),","), '[', 1))
dat.d7_8_neut$bin <- cut(dat.d7_8_neut$distance, breaks=seq(from=0,to=300, by=10))
dat.d7_8_neut$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d7_8_neut$bin),","), '[', 1))

dat.d1$bin <- cut(dat.d1$distance, breaks=seq(from=0,to=300, by=10))
dat.d1$bin1 <- gsub( '\\(', "", sapply(strsplit(as.character(dat.d1$bin),","), '[', 1))

a <- aggregate(mle_est ~ bin1 , data=dat.d7_7_sel, mean)
b <- aggregate(mle_est ~ bin1 , data=dat.d7_7_neut, mean)

c <- aggregate(mle_est ~ bin1 , data=dat.d7_8_sel, mean)
d <- aggregate(mle_est ~ bin1 , data=dat.d7_8_neut, mean)
e <- aggregate(mle_est ~ bin1 , data=dat.d1, mean)

f <- aggregate(mle_est ~ bin1 , data=dat.d7_d7_sel, mean)
g <- aggregate(mle_est ~ bin1 , data=dat.d7_d7_neut, mean)

png("~/urchin_af/figures/ld_binned_1.png", res=300, height=7, width=7, units="in")

plot(0,type='n', xlim=c(0,150), ylim=c(0,0.5),
    xlab="Distance between SNPs", ylab="Binned average estimated LD",
    main="Decay in LD with distance")

points(x=as.numeric(d$bin1), y=as.numeric(e$mle_est), col="black", pch=22, cex=1.5, lwd=1.5)

points(x=as.numeric(b$bin1), y=as.numeric(b$mle_est), col="red", cex=1.5, lwd=1.5)
points(x=as.numeric(d$bin1), y=as.numeric(d$mle_est), col="blue", pch=2, cex=1.5, lwd=1.5)

points(x=as.numeric(c$bin1), y=as.numeric(c$mle_est), col="blue", pch=17, cex=1.5)
points(x=as.numeric(a$bin1), y=as.numeric(a$mle_est), col="red", pch=19, cex=1.5)

points(x=as.numeric(g$bin1), y=as.numeric(g$mle_est), col="orange", pch=5, cex=1.5)
points(x=as.numeric(f$bin1), y=as.numeric(f$mle_est), col="orange", pch=18, cex=1.5)



legend(75,0.5,legend=c("pH day 7: selected","pH day 7: neutral", "Control day 7: selected", "Control day 7: neutral", "Control day 1: neutral"),
       col=c("red", "red", "blue", "blue", "black"), pch=c(19,1,17,2,0), cex=1.3)
dev.off()
