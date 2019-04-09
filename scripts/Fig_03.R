# figure 3

# a is LD decay, 
# b is starting af, folded

library(scales)
library(plyr)
library(dplyr)
#library(broom)
#library(tidyr)

setwd("~/urchin_af/analysis")
filelist = list.files(pattern = "^OASV2_DNA_*")
data_list = lapply(filelist, read.table, header=FALSE)

file_name <- gsub(".ld.out", "", filelist)
file_name <- gsub("OASV2_DNA_", "", file_name)

names(data_list)<-file_name
colnames <- c("chr", "snp1", "snp2", "x_11", "x_12", "x_21", "x_22", "af_a", "af_b", "dp_snp1", "dp_snp2","dp_intersect", "mle_low", "mle_est", "mle_hi", "r2", "a1", "a2", "b1", "b2")
dat <- lapply(data_list, setNames, colnames)

dat <- lapply(data_list, setNames, colnames)

dat <- lapply(dat, function(x) transform(x, id1 = paste(chr, snp1, sep=":")))
dat <- lapply(dat, function(x) transform(x, id2 = paste(chr, snp2, sep=":")))

#read in cmh results
cmh <-read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

cmh$sel_pH7 <- cmh$pH7_sig
cmh$sel_pH8<- cmh$pH8_sig

cmh$SNP <- paste(cmh$CHROM, cmh$POS, sep=":")

selected <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, cmh$sel_pH7, cmh$sel_pH8 )
colnames(selected) <- c("chr", "pos", "id", "sel_pH7", "sel_pH8")

df <- list()
for (i in 1:length(dat)){

    df[[i]] <- merge(dat[[i]], selected, by.x="id1", by.y="id")
    df[[i]] <- merge(df[[i]], selected, by.x="id2", by.y="id")    

}

##############
###
### plot decay in ld
###
##############
    
# calc distance between snp1 and snp2
df <- lapply(df, function(x) transform(x, distance = abs(x$snp1-x$snp2)))

df.1 <- df

df <- list()
for (i in 1:length(df.1)){

    df[[i]] <- df.1[[i]][which(df.1[[i]]$distance <=100),] 

}

names(df) <- names(dat)

# add read depth filter
df.1 <- df

df <- list()
for (i in 1:length(df.1)){

    df[[i]] <- df.1[[i]][which(df.1[[i]]$dp_intersect >10),] 

}

names(df) <- names(dat)

# d1 samples

d1 <- df[grep("D1_8_0", names(df))]

#subset to day 7 ph7 samples, pull out selected and neutral

d7_7 <- df[grep("D7_7_5", names(df))]

d7_7_sel <- list()
for (i in 1:length(d7_7)){

    d7_7_sel[[i]] <- d7_7[[i]][which(d7_7[[i]]$sel_pH7.x == TRUE | d7_7[[i]]$sel_pH7.y == TRUE ),]

}

d7_7_neut <- list()
for (i in 1:length(d7_7)){

    d7_7_neut[[i]] <- d7_7[[i]][which(d7_7[[i]]$sel_pH7.x == FALSE & d7_7[[i]]$sel_pH7.y == FALSE ),]

}

#subset to day 7 ph8 samples, pull out selected and neutral

d7_8 <- df[grep("D7_8_0", names(df))]

d7_8_sel <- list()
for (i in 1:length(d7_8)){

    d7_8_sel[[i]] <- d7_8[[i]][which(d7_8[[i]]$sel_pH8.x == TRUE | d7_8[[i]]$sel_pH8.y == TRUE ),]

}

d7_8_neut <- list()
for (i in 1:length(d7_8)){

    d7_8_neut[[i]] <- d7_8[[i]][which(d7_8[[i]]$sel_pH8.x == FALSE & d7_8[[i]]$sel_pH8.y == FALSE ),]

}

# merge all and perform stats

d1 <- lapply(d1, function(x) transform(x, pH = c("control")))
d7_7_sel <- lapply(d7_7_sel, function(x) transform(x, pH = c("pH_7.5")))
d7_7_neut <- lapply(d7_7_neut, function(x) transform(x, pH = c("pH_7.5")))
d7_8_sel <- lapply(d7_8_sel, function(x) transform(x, pH = c("pH_8.0")))
d7_8_neut <- lapply(d7_8_neut, function(x) transform(x, pH = c("pH_8.0")))

d1 <- lapply(d1, function(x) transform(x, group = c("neutral")))
d7_7_sel <- lapply(d7_7_sel, function(x) transform(x, group = c("selected")))
d7_7_neut <- lapply(d7_7_neut, function(x) transform(x, group = c("neutral")))
d7_8_sel <- lapply(d7_8_sel, function(x) transform(x, group = c("selected")))
d7_8_neut <- lapply(d7_8_neut, function(x) transform(x, group = c("neutral")))

d7_7_sel[[1]]$rep <- c("Rep_1")
d7_7_sel[[2]]$rep <- c("Rep_2")
d7_7_sel[[3]]$rep <- c("Rep_3")
d7_7_sel[[4]]$rep <- c("Rep_4")

d7_8_sel[[1]]$rep <- c("Rep_1")
d7_8_sel[[2]]$rep <- c("Rep_2")
d7_8_sel[[3]]$rep <- c("Rep_3")
d7_8_sel[[4]]$rep <- c("Rep_4")

d7_8_neut[[1]]$rep <- c("Rep_1")
d7_8_neut[[2]]$rep <- c("Rep_2")
d7_8_neut[[3]]$rep <- c("Rep_3")
d7_8_neut[[4]]$rep <- c("Rep_4")

d7_7_neut[[1]]$rep <- c("Rep_1")
d7_7_neut[[2]]$rep <- c("Rep_2")
d7_7_neut[[3]]$rep <- c("Rep_3")
d7_7_neut[[4]]$rep <- c("Rep_4")

d1[[1]]$rep <- c("Rep_1")
d1[[2]]$rep <- c("Rep_2")
d1[[3]]$rep <- c("Rep_3")
d1[[4]]$rep <- c("Rep_4")

d7_7_sel_all <- rbind.fill(d7_7_sel)
d7_8_sel_all <- rbind.fill(d7_8_sel)
d7_7_neut_all <- rbind.fill(d7_7_neut)
d7_8_neut_all <- rbind.fill(d7_8_neut)
d1_all <- rbind.fill(d1)

all_dat_out <- rbind(d7_7_sel_all, d7_8_sel_all)
all_dat_out <- rbind(all_dat_out, d7_7_neut_all)
all_dat_out <- rbind(all_dat_out, d7_8_neut_all)
all_dat_out <- rbind(all_dat_out, d1_all)

# difference in slopes
library(nlme)

m1 <- lme(mle_est ~ log10(distance),random=~1|rep,data=d7_7_sel_all)
anova(m1)
summary(m1)
# intercept: .31; distance: -0.06874804
m1 <- lme(mle_est ~ log10(distance),random=~1|rep,data=d7_8_sel_all)
anova(m1)
summary(m1)
# intercept:  0.2826; distance: -0.07159262
m1 <- lme(mle_est ~ log10(distance),random=~1|rep,data=d7_7_neut_all)
anova(m1)
summary(m1)
# intercept: .34; distance: -0.1127530
m1 <- lme(mle_est ~ log10(distance),random=~1|rep,data=d7_8_neut_all)
anova(m1)
summary(m1)
# intercept: .34; distance: -0.1126124
m1 <- lme(mle_est ~ log10(distance),random=~1|rep,data=d1_all)
anova(m1)
summary(m1)
# intercept: .34; distance: -0.1155936


## sample random subgroups of data. compare to selected loci

## sample random subgroups of data. compare to selected loci

perm_length <- length(which(cmh$sel_pH7 == TRUE))

perm.out <- list()
perm.raw <- list()

for (j in 1:4){
    for (i in 1:125){
        rand.sel <- rep("FALSE", nrow(cmh))
        rand.sel[sample(1:nrow(cmh), perm_length)]<- TRUE
        selected.perm <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, rand.sel )
        colnames(selected.perm) <- c("chr", "pos", "id", "selected")
        dat.perm <- merge(d7_7[[j]], selected.perm, by.x="id1", by.y="id", all.x=TRUE)
        dat.perm <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
        dat.perm$distance <- abs(dat.perm$snp1-dat.perm$snp2)
        dat.perm <- dat.perm[which(dat.perm$distance <=100),]
        dat.perm <- dat.perm[which(dat.perm$selected.y == TRUE | dat.perm$selected.x == TRUE),]
        mod.perm <- lm(mle_est ~ log10(distance), data = dat.perm)
        #lines(sort(dat.perm$distance, decreasing=FALSE),
        #    sort(fitted(mod.perm), decreasing=TRUE),
        #    lwd=1, col=rgb(0,0,0,alpha=0.1))
        if (i%%25 == 0){print(paste("current i loop",i))} # printing progress
        perm.out[[(length(perm.out)+1)]] <- cbind(sort(dat.perm$distance, decreasing=FALSE), sort(fitted(mod.perm), decreasing=TRUE))
        perm.raw[[(length(perm.out)+1)]] <- dat.perm
    
    }
    print(paste("Done with j loop:  ", j))
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
                            breaks=seq(from=0, to=100, by=1),
                            labels=seq(from=1, to=100, by=1)))
    out[[i]] <- data.frame(bin=seq(from=1, to=100, by=1), ld=tapply(tmp.df$y,tmp.df$bin , mean))
}

out.new <- do.call(rbind, out)

densities.qtiles <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(ld, 0.025, na.rm=TRUE),
            q50 = quantile(ld, 0.5, na.rm=TRUE),
            q95 = quantile(ld, 0.975, na.rm=TRUE))

#compare out.new to selected and non-selected data
## estimate intercept and slope of all perms

intercept <- c()
slope <- c()

for (i in 2:length(perm.raw)){
    m1 <- lm(mle_est ~ log10(distance),data=perm.raw[[i]])
    intercept[i] <- summary(m1)$coefficients[1,1]
    slope[i] <- summary(m1)$coefficients[2,1]

}

quantile(intercept, c(0.025, 0.975), na.rm=TRUE)
quantile(slope, c(0.025, 0.975), na.rm=TRUE)

#> quantile(intercept, c(0.025, 0.975), na.rm=TRUE)
#     2.5%     97.5%
#0.2977177 0.3816689
#> quantile(slope, c(0.025, 0.975), na.rm=TRUE)
#       2.5%       97.5%
#-0.13800580 -0.08461048


##################################################
##
## 8.0 perm
##
#################################

perm_length <- length(which(cmh$sel_pH8 == TRUE))

perm.out2 <- list()
perm.raw2 <- list()

for (j in 1:4){
    for (i in 1:125){
        rand.sel <- rep("FALSE", nrow(cmh))
        rand.sel[sample(1:nrow(cmh), perm_length)]<- TRUE
        selected.perm <- data.frame( cmh$CHROM, cmh$POS, cmh$SNP, rand.sel )
        colnames(selected.perm) <- c("chr", "pos", "id", "selected")
        dat.perm <- merge(d7_8[[j]], selected.perm, by.x="id1", by.y="id", all.x=TRUE)
        dat.perm <- merge(dat.perm, selected.perm, by.x="id2", by.y="id", all.x=TRUE)
        dat.perm$distance <- abs(dat.perm$snp1-dat.perm$snp2)
        dat.perm <- dat.perm[which(dat.perm$distance <=100),]
        dat.perm <- dat.perm[which(dat.perm$selected.y == TRUE | dat.perm$selected.x == TRUE),]
        mod.perm <- lm(mle_est ~ log10(distance), data = dat.perm)
        #lines(sort(dat.perm$distance, decreasing=FALSE),
        #    sort(fitted(mod.perm), decreasing=TRUE),
        #    lwd=1, col=rgb(0,0,0,alpha=0.1))
        if (i%%25 == 0){print(paste("current i loop",i))} # printing progress
        perm.out2[[(length(perm.out2)+1)]] <- cbind(sort(dat.perm$distance, decreasing=FALSE), sort(fitted(mod.perm), decreasing=TRUE))
        perm.raw2[[(length(perm.out2)+1)]] <- dat.perm
    
    }
    print(paste("Done with j loop:  ", j))
}



# for each permutation, take the mean of the bin, then combine all perms.

out <- list()
for (i in 1: length(perm.out2)){

    tmp.df <- data.frame(x=perm.out2[[i]][,1],
                            y= perm.out2[[i]][,2],
                        bin = cut(
                            perm.out2[[i]][,1],
                            breaks=seq(from=0, to=100, by=1),
                            labels=seq(from=1, to=100, by=1)))
    out[[i]] <- data.frame(bin=seq(from=1, to=100, by=1), ld=tapply(tmp.df$y,tmp.df$bin , mean))
}

out.new <- do.call(rbind, out)

densities.qtiles <- out.new %>%
group_by(bin) %>%
  summarise(q05 = quantile(ld, 0.025, na.rm=TRUE),
            q50 = quantile(ld, 0.5, na.rm=TRUE),
            q95 = quantile(ld, 0.975, na.rm=TRUE))

ld.qtiles <- densities.qtiles

#compare out.new to selected and non-selected data
## estimate intercept and slope of all perms

intercept <- c()
slope <- c()

for (i in 2:length(perm.raw2)){
    m1 <- lm(mle_est ~ log10(distance),data=perm.raw2[[i]])
    intercept[i] <- summary(m1)$coefficients[1,1]
    slope[i] <- summary(m1)$coefficients[2,1]

}

quantile(intercept, c(0.025, 0.975), na.rm=TRUE)
quantile(slope, c(0.025, 0.975), na.rm=TRUE)


# calc indiv models for plotting:

mod.d1  <- lapply(d1, function(x) lm(mle_est ~ log10(distance), data = x))
mod.d7_7_sel   <- lapply(d7_7_sel,  function(x) lm(mle_est ~ log10(distance), data = x))
mod.d7_7_neut  <- lapply(d7_7_neut, function(x) lm(mle_est ~ log10(distance), data = x))
mod.d7_8_sel   <- lapply(d7_8_sel,  function(x) lm(mle_est ~ log10(distance), data = x))
mod.d7_8_neut  <- lapply(d7_8_neut, function(x) lm(mle_est ~ log10(distance), data = x))


tiff("~/urchin_af/figures/Fig_03.tiff", height=85, width=170, units="mm", res=300)
par(mfrow = c(1, 2), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(0,type='n', xlim=c(1,100), ylim=c(0,0.38),
    main="",
    ylab="",
    xlab="",
    cex.lab=0.9, cex.axis=0.7,
    xaxt="n",yaxt="n")

polygon(x=c(densities.qtiles$bin,rev(densities.qtiles$bin)),
    y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
    col=alpha("black", alpha=0.3),border=NA)

lines(densities.qtiles$bin, densities.qtiles$q50, lwd=3, lty=2, col='black')


lines(sort(d7_7_neut[[1]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut[[1]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='firebrick3')
lines(sort(d7_7_neut[[2]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut[[2]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='firebrick3')
lines(sort(d7_7_neut[[3]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut[[3]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='firebrick3')
lines(sort(d7_7_neut[[4]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_neut[[4]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='firebrick3')

lines(sort(d7_8_neut[[1]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut[[1]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='royalblue3')
lines(sort(d7_8_neut[[2]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut[[2]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='royalblue3')
lines(sort(d7_8_neut[[3]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut[[3]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='royalblue3')
lines(sort(d7_8_neut[[4]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_neut[[4]]), decreasing=TRUE),
    lwd=1.5, lty=2, col='royalblue3')


lines(sort(d1[[1]]$distance, decreasing=FALSE), sort(fitted(mod.d1[[1]]), decreasing=TRUE),
    lwd=1.5, col='black', lty=2)
lines(sort(d1[[2]]$distance, decreasing=FALSE), sort(fitted(mod.d1[[2]]), decreasing=TRUE),
    lwd=1.5, col='black', lty=2)
lines(sort(d1[[3]]$distance, decreasing=FALSE), sort(fitted(mod.d1[[3]]), decreasing=TRUE),
    lwd=1.5, col='black', lty=2)
lines(sort(d1[[4]]$distance, decreasing=FALSE), sort(fitted(mod.d1[[4]]), decreasing=TRUE),
    lwd=1.5, col='black', lty=2)

# selected sites
lines(sort(d7_8_sel[[1]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel[[1]]), decreasing=TRUE),
    lwd=1.5, col='royalblue3')
lines(sort(d7_8_sel[[2]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel[[2]]), decreasing=TRUE),
    lwd=1.5, col='royalblue3')
lines(sort(d7_8_sel[[3]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel[[3]]), decreasing=TRUE),
    lwd=1.5, col='royalblue3')
lines(sort(d7_8_sel[[4]]$distance, decreasing=FALSE), sort(fitted(mod.d7_8_sel[[4]]), decreasing=TRUE),
    lwd=1.5, col='royalblue3')

lines(sort(d7_7_sel[[1]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel[[1]]), decreasing=TRUE),
    lwd=1.5, col='firebrick3')
lines(sort(d7_7_sel[[2]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel[[2]]), decreasing=TRUE),
    lwd=1.5, col='firebrick3')
lines(sort(d7_7_sel[[3]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel[[3]]), decreasing=TRUE),
    lwd=1.5, col='firebrick3')
lines(sort(d7_7_sel[[4]]$distance, decreasing=FALSE), sort(fitted(mod.d7_7_sel[[4]]), decreasing=TRUE),
    lwd=1.5, col='firebrick3')

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7, tcl=-0.2, at=c(1, 20, 40, 60, 80,100)) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Distance between SNPs in base pairs", line=1.5, cex.lab=0.9)
title(ylab="Estimated Linkage Disequilibrium", line=1.5, cex.lab=0.9)


legend("bottomleft",
    legend=c(expression('T'[0]),
            "pH 7.5: selected",
            "pH 7.5: neutral",
            "pH 8.0: selected",
            "pH 8.0: neutral"),
       col=c("black", "firebrick3", "firebrick3", "royalblue3", "royalblue3"),
       lty=c(2,1,2,1,2), lwd=2.2, cex=0.75, bty = "n")


mtext(text=bquote(paste('(',italic('a'),')')),
              side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)

#########
#
# unfolded maf
#
#########


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

#png("~/urchin_af/figures/Fig_S04_unfolded_af.png", height=85, width=85, units="mm", res=300)

plot(density(0:1), ylim=c(0,4),xlim=c(0,0.5), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(xlab="Starting allele frequency", line=1.5, cex.lab=0.9)
title(ylab="Density", line=1.5, cex.lab=0.9)
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)



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


mtext(text=bquote(paste('(',italic('b'),')')),
              side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)

dev.off()


