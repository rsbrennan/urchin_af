
# looking at variance in allele frequencies
# in selected vs non-selected. 
# wa
# looking at variance in allele frequencies
# in selected vs non-selected. 
# want to see if drift stronger in D7 control, drives lack of results

#read in cmh results
cmh <-read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

cut_off <- quantile(cmh$control_selection_pval, 0.01, na.rm=TRUE)

cut_off
cmh$selected <- FALSE
cmh$d7_sel<- FALSE
cmh$lab_sel <- FALSE

cmh$selected[which(cmh$pH_selection_pval < cut_off)] <- TRUE
cmh$d7_sel[which(cmh$d7_selection_pval < cut_off)] <- TRUE

sel_d7_7 <- cmh[which(cmh$d7_sel == TRUE),]
nuet_d7_7 <- cmh[which(cmh$d7_sel == FALSE),]

sel_7 <- cmh[which(cmh$selected == TRUE),]
neut_7 <- cmh[which(cmh$selected == FALSE),]

dat.7_sel <- sel_7[,grep("_af", colnames(sel_7))]
dat.7_sel <- dat.7_sel[,grep("D7_7", colnames(dat.7_sel))]

dat.7_neut <- neut_7[,grep("_af", colnames(neut_7))]
dat.7_neut <- dat.7_neut[,grep("D7_7", colnames(dat.7_neut))]

dat.8_sel <- sel_7[,grep("_af", colnames(sel_7))]
dat.8_sel <- dat.7_sel[,grep("D7_8", colnames(dat.8_sel))]

dat.8_neut <- neut_7[,grep("_af", colnames(neut_7))]
dat.8_neut <- dat.8_neut[,grep("D7_8", colnames(dat.8_neut))]

dat.d7_8 <- cmh[,grep("_af", colnames(cmh))]
dat.d7_8 <- dat.d7_8[,grep("D7_8", colnames(dat.d7_8))]

dat.d1_8 <- cmh[,grep("_af", colnames(cmh))]
dat.d1_8 <- dat.d1_8[,grep("D1_8", colnames(dat.d1_8))]

dat.d7_7 <- cmh[,grep("_af", colnames(cmh))]
dat.d7_7 <- dat.d7_7[,grep("D7_7", colnames(dat.d7_7))]

hist(apply(dat.7_sel ,1,sd), breaks=30, col=rgb(0,0,0,alpha=0.2), freq=FALSE)
hist(apply(dat.7_neut ,1,sd), breaks=50, col=rgb(1,0,0,alpha=0.2), freq=FALSE, add=T)
hist(apply(dat.8_sel, 1 ,sd), breaks=30, col=rgb(0,1,0,alpha=0.2), freq=FALSE, add=T)
hist(apply(dat.8_neut, 1 ,sd), breaks=50, col=rgb(0,0,1,alpha=0.2), freq=FALSE, add=T)

#make list of afs:
dat.8_sel.sd <- apply(dat.8_sel, 1 ,sd)
dat.8_neut.sd <- apply(dat.8_neut, 1 ,sd)
dat.7_sel.sd <- apply(dat.7_sel, 1 ,sd)
dat.7_neut.sd <- apply(dat.7_neut, 1 ,sd)
dat.d1.sd <- apply(dat.d1, 1 ,sd)
z <- c("dat.8_sel.sd", "dat.8_neut.sd","dat.7_sel.sd", "dat.7_neut.sd","dat.d1.sd")
dataList <- lapply(z, get, envir=environment())
names(dataList) <- z
boxplot(dataList)

ks.test

ks.test(dat.7_sel.sd, dat.8_sel.sd)

dat.d7_8.sd <- apply(dat.d7_8, 1 ,sd)
dat.d1_8.sd <- apply(dat.d1_8, 1 ,sd)
dat.d7_7.sd <- apply(dat.d7_7, 1 ,sd)
z <- c("dat.d7_8.sd", "dat.d1_8.sd","dat.d7_7.sd")
dataList <- lapply(z, get, envir=environment())
names(dataList) <- c("day 7 pH 8", "day 1 pH 8", "day 7 pH 7")

png("~/urchin_af/figures/af_sd_1.png", res=300, height=7, width=10, units="in")
boxplot(dataList, ylab="SD of allele frequency")
dev.off()

ks.test(dat.d7_8.sd, dat.d1_8.sd)
ks.test(dat.d7_8.sd, dat.d7_7.sd)
ks.test(dat.d7_7.sd, dat.d1_8.sd)


#### compare sig af in d7, d1
dat.8_d1_sel <- sel_7[,grep("_af", colnames(sel_7))]
dat.8_d1_sel <- dat.8_d1_sel[,grep("D1_8", colnames(dat.8_d1_sel))]
dat.8_d1_sel.sd <- apply(dat.8_d1_sel ,1,sd)

dat.7_d7_sel <- sel_7[,grep("_af", colnames(sel_7))]
dat.7_d7_sel <- dat.7_d7_sel[,grep("D7_7", colnames(dat.7_d7_sel))]
dat.7_d7_sel.sd <- apply(dat.7_d7_sel ,1,sd)

dat.8_d7_sel <- sel_7[,grep("_af", colnames(sel_7))]
dat.8_d7_sel <- dat.8_d7_sel[,grep("D7_8", colnames(dat.8_d7_sel))]
dat.8_d7_sel.sd <- apply(dat.8_d7_sel ,1,sd)

hist(dat.8_d1_sel.sd, breaks=30, col=rgb(0,0,0,alpha=0.2), freq=FALSE)
hist(dat.7_d7_sel.sd, breaks=50, col=rgb(1,0,0,alpha=0.2), freq=FALSE, add=T)
hist(dat.8_d7_sel.sd, breaks=30, col=rgb(0,1,0,alpha=0.2), freq=FALSE, add=T)

#boxplot
zx <- c("dat.8_d1_sel.sd", "dat.7_d7_sel.sd","dat.8_d7_sel.sd")
dataList <- lapply(zx, get, envir=environment())
names(dataList) <- zx
boxplot(dataList)

ks.test(dat.8_d7_sel.sd, dat.7_d7_sel.sd)
ks.test(dat.8_d7_sel.sd,dat.8_d1_sel.sd)
ks.test(dat.7_d7_sel.sd,dat.8_d1_sel.sd)

###########
##
## check that we're not just pulling out low sd variants
##
###########

a <- c()
b <- c()
c <- c()

d7_7_vs_d1_8 <- c()
d7_8_vs_d1_8 <- c()

for (i in 1:2000){

    rand.sel <- cmh[,grep("_af", colnames(cmh))]
    rand.sel <- rand.sel[sample(1:nrow(cmh), nrow(sel_7)),]
    rand.sel <- rand.sel[,grep("D1_8", colnames(rand.sel))]
    rand.sel.sd <- apply(rand.sel,1,sd)
    # save median
    a[i] <- median(rand.sel.sd)

    # ks test
    d7_7_vs_d1_8[i] <- ks.test(rand.sel.sd, dat.8_d7_sel.sd)$p.value
    d7_8_vs_d1_8[i] <- ks.test(rand.sel.sd, dat.7_d7_sel.sd)$p.value

    rand.sel <- cmh[,grep("_af", colnames(cmh))]
    rand.sel <- rand.sel[sample(1:nrow(cmh), nrow(sel_7)),]
    rand.sel <- rand.sel[,grep("D7_8", colnames(rand.sel))]
    rand.sel.sd <- apply(rand.sel,1,sd)
    # save median
    b[i] <- median(rand.sel.sd)

    rand.sel <- cmh[,grep("_af", colnames(cmh))]
    rand.sel <- rand.sel[sample(1:nrow(cmh), nrow(sel_7)),]
    rand.sel <- rand.sel[,grep("D7_7", colnames(rand.sel))]
    rand.sel.sd <- apply(rand.sel,1,sd)
    # save median
    c[i] <- median(rand.sel.sd) 

    if (i%%100 == 0){print(i)} # printing progress

}


hist(d7_7_vs_d1_8, breaks=30, col=rgb(0,0,0,alpha=0.2), freq=FALSE)

png("~/urchin_af/figures/af_sd_2.png", res=300, height=7, width=10, units="in")
hist(a, breaks=40, col=rgb(0,0,0,alpha=0.2), freq=FALSE, xlim=c(0.032, 0.055), 
    main="Permuted median allele freq of 818 SNPs\nlines show median AF of 818 sig SNPs in each set of samples")
hist(b, breaks=40, col=rgb(0,1,0,alpha=0.2), freq=FALSE, xlim=c(0.032, 0.055), add=T)
hist(c, breaks=40, col=rgb(1,0,0,alpha=0.2), freq=FALSE, xlim=c(0.032, 0.055), add=T)

abline(v=median(dat.8_d1_sel.sd), col="grey", lwd=2)
abline(v=median(dat.8_d7_sel.sd), col="green", lwd=2)
abline(v=median(dat.7_d7_sel.sd), col="red", lwd=2)

legend(0.038, 500, legend=c("Day 1 pH 8","Day 7 pH 8", "Day 7 pH 7.5"),
       col=c("grey", "green", "red"), pch=c(19,19,19), cex=1.3)

dev.off()
