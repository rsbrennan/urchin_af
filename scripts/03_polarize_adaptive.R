
################################################################
################################################################
## script to polarize allele frequency based on "adaptive allele"
## then look at the starting af of the adaptive allele.
## Basically, whichever is increasing in response to low pH is adaptive.
################################################################
################################################################

library(dplyr)
library(reshape)
library(ggplot2)
library(qvalue)
library(scales)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

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
dat_snp <- paste(mydata$CHROM, mydata$POS, sep=":")
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

# next, want to figure out which allele is increasing in frequency from D1_8 to D7

# rm depth columns
af.1 <-  af[,grep("_DPtotal", colnames(af), invert=TRUE)]

# calculate mean of each group
#### calc AF for each rep

# calculate mean of each group
gp <- c("D1_8", "D7_7", "D7_8")
af.mean <-  data.frame(matrix(ncol=6, nrow=nrow(af.1)))
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

# figure out which are "adaptive". THat is, sig for either pH7 or pH8 or both. 
#If neither, this is arbitrary assignment. use mean from both D7


# next, want to figure out which allele is increasing in frequency from D1_8 to D7_7

af1_7 <- af.mean$D7_7_DP1 - af.mean$D1_8_DP1
af2_7 <- af.mean$D7_7_DP2 - af.mean$D1_8_DP2
af1_8 <- af.mean$D7_8_DP1 - af.mean$D1_8_DP1
af2_8 <- af.mean$D7_8_DP2 - af.mean$D1_8_DP2
af1_both <- rowMeans(cbind(af1_7, af1_8))

af_out <- as.data.frame(matrix(nrow=nrow(af.mean), ncol=3))
colnames(af_out) <- c("D1_8_af", "D7_7_af", "D7_8_af") 
# mydata is in same order as af.mean. so can check sig cols. 
# let's also calculate the change in allele freq based on each rep, rather than mean.


af_out_in <- as.data.frame(matrix(nrow=nrow(af.mean), ncol=16))
colnames(af_out_in) <- colnames(af.1[9:ncol(af.1)])

af_out_reps <- as.data.frame(matrix(nrow=nrow(af.mean), ncol=8))
colnames(af_out_reps) <- unique(substr(colnames(af_out_in), 1,7))

af_out_in$D7_7_19_DP1 <- af.1$D7_7_19_DP1 - af.1$D1_8_01_DP1 
af_out_in$D7_7_19_DP2 <- af.1$D7_7_19_DP2 - af.1$D1_8_01_DP2 
af_out_in$D7_7_21_DP1 <- af.1$D7_7_21_DP1 - af.1$D1_8_02_DP1 
af_out_in$D7_7_21_DP2 <- af.1$D7_7_21_DP2 - af.1$D1_8_02_DP2 
af_out_in$D7_7_22_DP1 <- af.1$D7_7_22_DP1 - af.1$D1_8_04_DP1 
af_out_in$D7_7_22_DP2 <- af.1$D7_7_22_DP2 - af.1$D1_8_04_DP2 
af_out_in$D7_7_23_DP1 <- af.1$D7_7_23_DP1 - af.1$D1_8_07_DP1 
af_out_in$D7_7_23_DP2 <- af.1$D7_7_23_DP2 - af.1$D1_8_07_DP2

af_out_in$D7_8_14_DP1 <- af.1$D7_8_14_DP1 - af.1$D1_8_01_DP1 
af_out_in$D7_8_14_DP2 <- af.1$D7_8_14_DP2 - af.1$D1_8_01_DP2 
af_out_in$D7_8_15_DP1 <- af.1$D7_8_15_DP1 - af.1$D1_8_02_DP1 
af_out_in$D7_8_15_DP2 <- af.1$D7_8_15_DP2 - af.1$D1_8_02_DP2 
af_out_in$D7_8_16_DP1 <- af.1$D7_8_16_DP1 - af.1$D1_8_04_DP1 
af_out_in$D7_8_16_DP2 <- af.1$D7_8_16_DP2 - af.1$D1_8_04_DP2 
af_out_in$D7_8_17_DP1 <- af.1$D7_8_17_DP1 - af.1$D1_8_07_DP1 
af_out_in$D7_8_17_DP2 <- af.1$D7_8_17_DP2 - af.1$D1_8_07_DP2



for (i in 1:nrow(af_out)) {
    if(mydata$pH7_sig[i] == TRUE & mydata$pH8_sig[i] == FALSE){
        if(af1_7[i] > 0){
            af_out[i,] <- af.mean[i, grep("DP1", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP1", colnames(af_out_in))]
        }
        else{
            af_out[i,] <- af.mean[i, grep("DP2", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP2", colnames(af_out_in))]

        }

    }
    if(mydata$pH8_sig[i] == TRUE & mydata$pH7_sig[i] == FALSE){
        if(af1_8[i] > 0){
            af_out[i,] <- af.mean[i, grep("DP1", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP1", colnames(af_out_in))]

        }
        else{
            af_out[i,] <- af.mean[i, grep("DP2", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP2", colnames(af_out_in))]

        }
    }
    if(mydata$pH8_sig[i] == FALSE & mydata$pH7_sig[i] == FALSE){
        if(af1_both[i] > 0){
            af_out[i,] <- af.mean[i, grep("DP1", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP1", colnames(af_out_in))]

        }
        else{
            af_out[i,] <- af.mean[i, grep("DP2", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP2", colnames(af_out_in))]

        }
    }
    if(mydata$pH8_sig[i] == TRUE & mydata$pH7_sig[i] == TRUE){
        if(af1_both[i] > 0){
            af_out[i,] <- af.mean[i, grep("DP1", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP1", colnames(af_out_in))]

        }
        else{
            af_out[i,] <- af.mean[i, grep("DP2", colnames(af.mean))]
            af_out_reps[i,] <- af_out_in[i, grep("DP2", colnames(af_out_in))]

        }
    }
}

colnames(af_out_reps) <- paste("delta",colnames(af_out_reps), sep="_")

delta_75 <- apply(af_out_reps[,1:4],1,mean)
delta_80 <- apply(af_out_reps[,5:8],1,mean)

sd_75 <- apply(af_out_reps[,1:4],1,sd)
sd_80 <- apply(af_out_reps[,5:8],1,sd)

# combine all:

out1 <- cbind(mydata, af_out)
out2 <- cbind(out1, af_out_reps)

out3 <- out2[ , !(names(out2) %in% c("D1_8_mean", "D7_7_mean", "D7_8_mean"))]

out3$mean_delta_75 <- delta_75
out3$mean_delta_80 <- delta_80
out3$sd_delta_75 <- sd_75
out3$sd_delta_80 <- sd_80

mean(out3$mean_delta_75[which(out3$pH7_sig == TRUE)])
mean(out3$sd_delta_75[which(out3$pH7_sig == TRUE)])
mean(out3$mean_delta_80[which(out3$pH8_sig == TRUE)])
mean(out3$sd_delta_80[which(out3$pH8_sig == TRUE)])
mean(out3$mean_delta_75[which(out3$pH7_sig == TRUE & out3$pH8_sig == TRUE)])
mean(out3$sd_delta_75[which(out3$pH7_sig == TRUE & out3$pH8_sig == TRUE)])

mean(out3$mean_delta_75)
mean(out3$sd_delta_75)

mean(out3$mean_delta_80)
mean(out3$sd_delta_75)


# save adaptive allele table

write.table(file="~/urchin_af/analysis/adaptive_allelefreq.txt",out3,
    col.name=TRUE, quote=FALSE, row.name=FALSE, sep= "\t" )
