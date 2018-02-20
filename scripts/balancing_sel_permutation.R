############################################
##
## look at starting af
## balancing selection
## permutation analysis
##
############################################

library(dplyr)
library(reshape)
library(ggplot2)
library(qvalue)
library(scales)

cmh_dat <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

snp.final <- af_out[which(cmh_dat$pH_sig == TRUE)]

# compare pH selection vs control


################################################################
################################################################
## script to polarize allele frequency based on "adaptive allele"
## then look at the starting af of the adaptive allele.
## Basically, whichever is increasing in response to low pH is adaptive.
################################################################
################################################################

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

snp.info <- read.table("~/urchin_af/analysis/af.info.txt", header=TRUE)

snp <- read.table(text= system("zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^##' | cut -f 10-", intern=TRUE),
    header=TRUE)
snp.info <- read.table(text= system("zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^##' | cut -f 1-8", intern=TRUE),
    header=FALSE)
colnames(snp.info) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

colnames(snp) <- gsub("OASV2_DNA_", "", colnames(snp))
colnames(snp) <- gsub("S_", "", colnames(snp))
colnames(snp) <- gsub("5_", "", colnames(snp))
colnames(snp) <- gsub("0_", "", colnames(snp))

# first order both dataframe
# Sort by column index [1] then [3]
snp.1 <- snp[
  order( snp.info[,1], snp.info[,2] ),
  ]

mydata.1 <- mydata[
  order( mydata[,1], mydata[,2] ),
  ]

snp <- snp.1
mydata <- mydata.1

out.split <- apply(snp, 2, function(x) strsplit(x, ":"))

dat <- data.frame(row.names=seq(from=1, to=nrow(snp),by= 1))

# get depths
for(i in 1:length(names(out.split))){
    ct <- matrix(unlist(out.split[[i]]), ncol=4, byrow=TRUE)
    #
    dat[,paste(names(out.split)[i], "DPtotal", sep="_")] <- sapply(strsplit(ct[,3], ","), "[", 1)
    dat[,paste(names(out.split)[i], "DP1", sep="_")] <- sapply(strsplit(ct[,4], ","), "[", 1)
    dat[,paste(names(out.split)[i], "DP2", sep="_")] <- sapply(strsplit(ct[,4], ","), "[", 2)
}

#dep <-  dat[,grep("_DPtotal", colnames(dat), invert=TRUE)]

#convert all to numeric
dep <- data.frame(sapply(dat, function(x) as.numeric(as.character(x))))

# want to figure out which allele is increasing in frequency from D1_8 to D7_7

#first, calculate AF
# make empty data frame to hold af
af <- data.frame(matrix(ncol=ncol(dep), nrow=nrow(dep)))
pop <- unique(substr(colnames(dep),1,7))
colnames(af) <- colnames(dep)

for(i in 1:ncol(af)){
    # pull out DP1
    dp1 <- as.numeric(dep[,grep(paste(pop[i], "_DP1", sep=""), colnames(dep))])
    # DP2
    dp2 <- as.numeric(dep[,grep(paste(pop[i], "_DP2", sep=""), colnames(dep))])
    # total depth
    dptot <- as.numeric(dep[,grep(paste(pop[i], "_DPtotal", sep=""), colnames(dep))])
    #calc af
    dp1.freq <- dp1/dptot
    dp2.freq <- dp2/dptot
    # add to af
    af[,grep(paste(pop[i], "_DP1", sep=""), colnames(af))] <- dp1.freq
    af[,grep(paste(pop[i], "_DP2", sep=""), colnames(af))] <- dp2.freq
    af[,grep(paste(pop[i], "_DPtotal", sep=""), colnames(af))] <- dptot
}

# next, want to figure out which allele is increasing in frequency from D1_8 to D7_7

# rm depth columns
af.1 <-  af[,grep("_DPtotal", colnames(af), invert=TRUE)]


# calculate mean of each group
#### calc AF for each rep


# calculate mean of each group
gp <- c("D1_7", "D1_8", "D7_7", "D7_8")
af.mean <-  data.frame(matrix(ncol=8, nrow=nrow(af.1)))
nm1 <- paste(gp, "_DP1", sep="")
nm2 <- paste(gp, "_DP2", sep="")
nm3 <- c(nm1, nm2)
colnames(af.mean) <- nm3

for (i in 1:length(gp)){
    sub.gp <- af.1[,grep(gp[i], colnames(af.1))]

    dp1 <- sub.gp[, grep("_DP1", colnames(sub.gp))]
    dp1.mean <- apply(dp1,1,mean)
    dp2 <- sub.gp[, grep("_DP2", colnames(sub.gp))]
    dp2.mean <- apply(dp2,1,mean)

    # add means to af.mean
    af.mean[[paste(gp[i], "_DP1", sep="")]] <- dp1.mean
    af.mean[[paste(gp[i], "_DP2", sep="")]] <- dp2.mean
}

# next, want to figure out which allele is increasing in frequency from D1_8 to D7_7

af1 <- af.mean$D7_7_DP1 - af.mean$D1_8_DP1
af2 <- af.mean$D7_7_DP2 - af.mean$D1_8_DP2

af_out <- rep(NA, nrow(af.mean))
af_d7 <- rep(NA, nrow(af.mean))

for (i in 1:nrow(af.1)) {
    if(af1[i] > 0 ){
        af_out[i] <- af.mean$D1_8_DP1[i]
        af_d7[i] <- af.mean$D7_7_DP1[i]
    }
    else{
        af_out[i] <- af.mean$D1_8_DP2[i]
        af_d7[i] <- af.mean$D7_7_DP2[i]

    }
}

# need to pull out only selected alleles

snp.sel <- af_out[which(mydata$pH_sig == TRUE)]
snp.af_d7 <- af_d7[which(mydata$pH_sig == TRUE)]
snp.ctr <- af_out[which(mydata$control_selection_pval < cut_off)]

png("~/urchin_af/figures/maf_unfolded_hist.png", res=300, height=7, width=7, units="in")
par(mfrow = c(1, 1))
hist(af_out, freq=FALSE, ylim=c(0,4.6), col=alpha("black", alpha = 0.6), breaks=50,
    main="Unfolded MAF", xlab="allele frequency")
hist(snp.ctr, freq=FALSE, add=T, col=alpha("blue", alpha = 0.4), breaks=50)
hist(snp.sel, freq=FALSE, add=T, col=alpha("red", alpha = 0.4), breaks=50)
legend("topright", c("D1 pH8- all", "D7 pH 8-selected", "D7 pH 7-selected"), pch=19, col=c("black", "blue", "red"), )
dev.off()

snp.final <- snp.sel



################################################
######
## permute to pull out false positives and compare sfs
######
################################################

# permute samples. pull out "responsive" loci. calc summary stats
# try to run in parallel

library(foreach)
library(doParallel)
library(scales)
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-6) #not to overload your computer. This is setting up 10 different computing environments
registerDoParallel(cl)
clusterCall(cl, function() library(scales)) # need to load library for each node.


#png("~/urchin_af/figures/maf_unfolded_permute.png", res=301, height=7, width=7, units="in")

#plot(density(af_out), ylim=c(0,3), lwd=0,
 #   main="", xlab="allele frequency")

#ks.out <- c()

# this transposes the results
comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
}

my_results_par <- foreach(perm_rep = 1:1000, .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar%
    {
#    sink("~/monitor.txt",append=TRUE)
    cat(paste("Starting iteration",perm_rep,"\n"),
       file="~/log.monitor.txt", append=TRUE)

   control_selection_pval <- c()

    DP1 <- mydata[,grep("_DP1", colnames(mydata))]
    DP2 <- mydata[,grep("_DP2", colnames(mydata))]
    shuf <- sample(seq(1:ncol(DP1)))
    DP1 <- DP1[,shuf]
    DP2 <- DP2[,shuf]

    dp1.sub <- DP1[,1:8]
    colnames(dp1.sub) <- c("D1_1_dp1","D1_2_dp1","D1_3_dp1","D1_4_dp1","D7_1_dp1","D7_2_dp1","D7_3_dp1","D7_4_dp1")
    dp2.sub <- DP2[,1:8]
    colnames(dp2.sub) <- c("D1_1_dp2","D1_2_dp2","D1_3_dp2","D1_4_dp2","D7_1_dp2","D7_2_dp2","D7_3_dp2","D7_4_dp2")
    for(i in 1:nrow(dp2.sub)){

        sub_dp1 <- stack(dp1.sub[i,])
        sub_dp1$allele <- rep("ac1", nrow(sub_dp1))
        sub_dp2 <- stack(dp2.sub[i,])
        sub_dp2$allele <- rep("ac2", nrow(sub_dp2))
        colnames(sub_dp2) <- c("count", "ind", "allele")
        colnames(sub_dp1) <- c("count", "ind", "allele")
        sub_all <- rbind(sub_dp1, sub_dp2)
        sub_all$day <- substr(sub_all$ind, 1,2)
        sub_all$replicate <- substr(sub_all$ind, 4,4)
        Data.xtabs = xtabs(count ~ allele + day + replicate,
                       data=sub_all)
        test <- mantelhaen.test(Data.xtabs)
        control_selection_pval[i] <- test$p.value

        #if (i%%1000 == 0){print(i)}

        #ftable(Data.xtabs)
    }

    # pull out sfs of sig results

cut_off <- quantile(control_selection_pval, 0.001, na.rm=TRUE)
out <- cbind(dp1.sub,dp2.sub )
pop <- colnames(out)

af <- data.frame(matrix(ncol=ncol(out), nrow=nrow(out)))
colnames(af) <- pop

for(i in 1:8){
    #pull out pop lab
    pop_nm <- substr(pop[i],1,4)
    # pull out DP1
    dp1 <- as.numeric(out[,grep(paste(pop_nm, "_dp1", sep=""), colnames(out))])
    # total depth
    dp2 <- as.numeric(out[,grep(paste(pop_nm, "_dp2", sep=""), colnames(out))])
    #calc af
    a.freq <- dp1/(dp1+dp2)
    # add to af
    af[,i] <- a.freq
    #print(i)
}
for(i in 9:16){
    #pull out pop lab
    pop_nm <- substr(pop[i],1,4)
    # pull out DP1
    dp1 <- as.numeric(out[,grep(paste(pop_nm, "_dp1", sep=""), colnames(out))])
    # total depth
    dp2 <- as.numeric(out[,grep(paste(pop_nm, "_dp2", sep=""), colnames(out))])
    #calc af
    a.freq <- dp2/(dp1+dp2)
    # add to af
    af[,i] <- a.freq
    #print(i)
}

# calculate mean of each group
#### calc AF for each rep

# calculate mean of each group
gp <- c("D1", "D7")
af.mean <-  data.frame(matrix(ncol=4, nrow=nrow(af)))
nm1 <- paste(gp, "_dp1", sep="")
nm2 <- paste(gp, "_dp2", sep="")
nm3 <- c(nm1, nm2)
colnames(af.mean) <- nm3

for (i in 1:length(gp)){
    sub.gp <- af[,grep(gp[i], colnames(af))]

    dp1 <- sub.gp[, grep("_dp1", colnames(sub.gp))]
    dp1.mean <- apply(dp1,1,mean)
    dp2 <- sub.gp[, grep("_dp2", colnames(sub.gp))]
    dp2.mean <- apply(dp2,1,mean)

    # add means to af.mean
    af.mean[[paste(gp[i], "_dp1", sep="")]] <- dp1.mean
    af.mean[[paste(gp[i], "_dp2", sep="")]] <- dp2.mean
}

# next, want to figure out which allele is increasing in frequency from D1_8 to D7_7

af1 <- af.mean$D7_dp1 - af.mean$D1_dp1
af2 <- af.mean$D7_dp2 - af.mean$D1_dp2

af_out <- rep(NA, nrow(af.mean))
af_d7 <- rep(NA, nrow(af.mean))

for (i in 1:nrow(af)) {
    if(af1[i] > 0 ){
        af_out[i] <- af.mean$D1_dp1[i]
        af_d7[i] <- af.mean$D7_dp1[i]
    }
    else{
        af_out[i] <- af.mean$D1_dp2[i]
        af_d7[i] <- af.mean$D7_dp2[i]

    }
}

cut_off <- quantile(control_selection_pval, 0.01, na.rm=TRUE)

# need to pull out only selected alleles
snp.sel <- af_out[which(control_selection_pval < cut_off)]

#lines(density(snp.sel), col=alpha("black", 0.15), lwd=1)

# I think this might just add it to my_results_par

list(ks.test(snp.final,snp.sel)$p.value, snp.sel)

}


write.table(my_results_par[[1]], file="~/urchin_af/analysis/permutation_pval.txt",  col.names=TRUE, quote=FALSE, sep="\t")
write.table(my_results_par[[2]], file="~/urchin_af/analysis/permutation_af.txt",  col.names=TRUE, quote=FALSE, sep="\t")

#stop cluster
stopCluster(cl)

# 17/1000 > 0.05. 


png("~/urchin_af/figures/maf_unfolded_permute.png", res=301, height=7, width=7, units="in")

plot(density(af_out), ylim=c(0,3), lwd=0,
    main="", xlab="allele frequency")

for (i in 1: ncol(my_results_par[[2]])){
    lines(density(unlist(my_results_par[[2]][,i])), col=alpha("black", 0.1), lwd=1)
}

lines(density(snp.final), col=alpha("red", 1), lwd=3)

legend("topright", c("permuted", "D7 pH 7.5-selected"), pch=19, col=c("black", "red"))

dev.off()

# look at proportion significant

which(unlist(my_results_par[[1]]) < 0.05)

## plot results, calculate p values

png("~/urchin_af/figures/all_permutations.png", res=300, height=14, width=10, units="in")
par(mfrow = c(3, 2))

hist(my_results_par$perm.mean, col="grey", breaks=100, main="mean af",
    xlim=c(0.165, 0.175))
abline(v=mean(apply(sel.sfs,1,mean)), col="red", lwd=3)

hist(my_results_par$perm.median, col="grey", breaks=100, , main="median af",
    xlim=c(min(my_results_par$perm.median), .16))
abline(v=mean(apply(sel.sfs,1,median), 0.95), col="red", lwd=3)

hist(my_results_par$perm.75, col="grey", breaks=100, main="75 percentile",
    xlim=c(min(my_results_par$perm.75), 0.28))
abline(v=quantile(apply(sel.sfs,1,mean), 0.75), col="red", lwd=3)

hist(my_results_par$perm.90, col="grey", breaks=100, main="90 percentile",
    xlim=c(min(my_results_par$perm.90), quantile(apply(sel.sfs,1,mean), 0.90)*1.05))
abline(v=quantile(apply(sel.sfs,1,mean), 0.90), col="red", lwd=3)

hist(my_results_par$perm.95, col="grey", breaks=100, main="95 percentile",
    xlim=c(0.37, 0.45))
abline(v=quantile(apply(sel.sfs,1,mean), 0.95), col="red", lwd=3)

dev.off()
