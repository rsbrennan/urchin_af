
## 06_balancing_sel_permutation.R

# running the permutation for allele freq distributions

library(dplyr)
library(reshape)
library(ggplot2)
library(qvalue)
library(scales)

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", stringsAsFactors=FALSE, header=TRUE)

cut_off <- 0.05/(9828)
# need to pull out only selected alleles
snp.sel_75 <- mydata$D1_8_af[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off)]
snp.sel_80 <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off)]
snp.sel_both <- mydata$D1_8_af[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off)]

png("~/urchin_af/figures/maf_unfolded_hist.png", res=300, height=7, width=7, units="in")
par(mfrow = c(1, 1))
hist(mydata$D1_8_af, freq=FALSE, ylim=c(0,6), col=alpha("black", alpha = 0.6), breaks=50,
    main="Unfolded MAF", xlab="allele frequency")
hist(snp.sel_80, freq=FALSE, add=T, col=alpha("royalblue3", alpha = 0.4), breaks=50)
hist(snp.sel_75, freq=FALSE, add=T, col=alpha("firebrick3", alpha = 0.4), breaks=50)
hist(snp.sel_both, freq=FALSE, add=T, col=alpha("darkorchid4", alpha = 0.4), breaks=50)
legend("topright", c("D1 pH8- all", "D7 pH 8-selected", "D7 pH 7-selected"), pch=19, col=c("black", "blue", "red"), )
dev.off()


################################################
######
## permute to pull out false positives and compare sfs
######
################################################

# permute samples. pull out "responsive" loci. calc summary stats
# basically, randomly assigning samples to groups. calculating cmh, looking at af, etc.

# try to run in parallel

library(foreach)
library(doParallel)
library(scales)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(8) #not to overload your computer. This is setting up 10 different computing environments
registerDoParallel(cl)
clusterCall(cl, function() library(scales)) # need to load library for each node.
clusterCall(cl, function() library(qvalue)) # need to load library for each node.

# this transposes the results
comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
}

my_results_par <- foreach(perm_rep = 1:500, .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar%
    {
#    sink("~/monitor.txt",append=TRUE)
    cat(paste("Starting iteration",perm_rep,"\n"),
       file="~/log.balancingsel.txt", append=TRUE)

   control_selection_pval <- c()

    DP1 <- mydata[,grep("_DP1", colnames(mydata))]
    DP2 <- mydata[,grep("_DP2", colnames(mydata))]
    # mix up these samples
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

    }

# pull out sfs of sig results

out <- cbind(dp1.sub,dp2.sub )
pop <- colnames(out)

af <- data.frame(matrix(ncol=ncol(out), nrow=nrow(out)))
colnames(af) <- pop

# calc af
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

for (i in 1:nrow(af)) {
    if(af1[i] > 0 ){
        af_out[i] <- af.mean$D1_dp1[i]
    }
    else{
        af_out[i] <- af.mean$D1_dp2[i]

    }
}

control_selection_pval[which(is.na(control_selection_pval))] <- 1 # bc some are invariant
cut_off <- quantile(control_selection_pval, 0.01)

# need to pull out only selected alleles, take top 0.015 quantile
snp.sel <- af_out[which(control_selection_pval < cut_off)]

# the following will add to output
list(ks.test(snp.sel_80, snp.sel, alternative="greater")$p.value,
    ks.test(snp.sel_75,snp.sel, alternative="greater")$p.value,
    ks.test(snp.sel_both,snp.sel, alternative="greater")$p.value,
    snp.sel,
)

}

write.table(my_results_par[[2]], file="~/urchin_af/analysis/permutation_75_pval.txt",  col.names=TRUE, quote=FALSE, sep="\t")
write.table(my_results_par[[1]], file="~/urchin_af/analysis/permutation_8_pval.txt",  col.names=TRUE, quote=FALSE, sep="\t")
write.table(my_results_par[[3]], file="~/urchin_af/analysis/permutation_both_pval.txt",  col.names=TRUE, quote=FALSE, sep="\t")
write.table(my_results_par[[4]], file="~/urchin_af/analysis/permutation_af.txt",  col.names=TRUE, quote=FALSE, sep="\t")



#stop cluster
stopCluster(cl)

# which comparisons from the ks test are sig different?
length(which(my_results_par[[1]] < 0.05)) # ks, 8.0 vs perm
# 491
length(which(my_results_par[[2]] < 0.05)) # ks, 7.5 vs perm
#493
length(which(my_results_par[[3]] < 0.05)) # ks, both vs perm
# 412

png("~/urchin_af/figures/maf_unfolded_permute.png", res=300, height=7, width=7, units="in")

plot(density(0:1), ylim=c(0,3),xlim=c(0,1), lwd=0,
    main="", xlab="allele frequency")
avg_perm <- c()
for (i in 1: ncol(my_results_par[[4]])){
    lines(density(unlist(my_results_par[[4]][,i]),na.rm=TRUE), col=alpha("black", 0.08), lwd=1)
    avg_perm <- c(avg_perm, unlist(my_results_par[[4]][,i]))
}

lines(density(snp.sel_75), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both), col=alpha("darkorchid2", 1), lwd=3)

abline(v=mean(avg_perm), lty=2, col= "black", lwd=3)
abline(v=mean(snp.sel_75), lty=2, col= "firebrick3", lwd=3)
abline(v=mean(snp.sel_80), lty=2, col= "royalblue3", lwd=3)
abline(v=mean(snp.sel_both), lty=2, col= "darkorchid2", lwd=3)

legend("topright", c("permuted", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    pch=19, col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=NULL)

dev.off()

