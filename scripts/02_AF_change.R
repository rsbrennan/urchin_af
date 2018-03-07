library(dplyr)
library(reshape)
library(ggplot2)
library(qvalue)
library(scales)
library(qqman)

mydata <- read.table("~/urchin_af/data/allele.freq.txt", header=TRUE)

## subset data to only include af measures
dat <- mydata[,grep("_af", colnames(mydata))]

####################
##
## Cochran–Mantel–Haenszel Test
##
####################

# cmh for pH selection after 7 days, only include d1 8 ph as control
pH_selection_pval <-c()

for(i in 1:nrow(mydata)){

    sub_ad <- mydata[i,grep("_DP1", colnames(mydata))]
    sub_ad <- stack(sub_ad)
    sub_ad$allele <- rep("ac1", nrow(sub_ad))
    sub_dp <- mydata[i,grep("_DP2", colnames(mydata))]
    sub_dp <- stack(sub_dp)
    sub_dp$allele <- rep("ac2", nrow(sub_dp))
    colnames(sub_dp) <- c("count", "ind", "allele")
    colnames(sub_ad) <- c("count", "ind", "allele")
#   count <- sub_dp$depth - sub_ad$count
#   sub_ac2 <- data.frame(count=count, ind=sub_dp$ind, allele=sub_dp$allele)
    sub_all <- rbind(sub_ad, sub_dp)
    sub_all$day <- substr(sub_all$ind, 1,2)
    sub_all$replicate <- substr(sub_all$ind, 6,7)
    sub_all$pH <- substr(sub_all$ind, 4,4)

    # remove unwanted replicates
    sub_all <- sub_all[grep("D1_7_", sub_all$ind, invert=TRUE),]
    sub_all <- sub_all[grep("D7_8_", sub_all$ind, invert=TRUE),]
    # need to make replicates match. pairings are arbitrary at this point
    sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
            paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)

    Data.xtabs = xtabs(count ~ allele + day + replicate,
                   data=sub_all)
    test <- mantelhaen.test(Data.xtabs)
    pH_selection_pval[i] <- test$p.value

    if (i%%5000 == 0){print(i)}
    #ftable(Data.xtabs)
}

pH_selection_pval[which(is.na(pH_selection_pval))] <- 1 # bc some are invariant

print("D1-D7 low pH cmh done")

#############
# cmh for pH selection after 1 day1,  include d1 8 ph as control

d1_selection_pval <-c()

for(i in 1:nrow(mydata)){

    sub_ad <- mydata[i,grep("_DP1", colnames(mydata))]
    sub_ad <- stack(sub_ad)
    sub_ad$allele <- rep("ac1", nrow(sub_ad))
    sub_dp <- mydata[i,grep("_DP2", colnames(mydata))]
    sub_dp <- stack(sub_dp)
    sub_dp$allele <- rep("ac2", nrow(sub_dp))
    colnames(sub_dp) <- c("count", "ind", "allele")
    colnames(sub_ad) <- c("count", "ind", "allele")
    sub_all <- rbind(sub_ad, sub_dp)
    sub_all$day <- substr(sub_all$ind, 1,2)
    sub_all$replicate <- substr(sub_all$ind, 6,7)
    sub_all$pH <- substr(sub_all$ind, 4,4)
    # remove unwanted replicates
    sub_all <- sub_all[grep("D7_7_", sub_all$ind, invert=TRUE),]
    sub_all <- sub_all[grep("D7_8_", sub_all$ind, invert=TRUE),]
    sub_all <- sub_all[grep("D1_8_07", sub_all$ind, invert=TRUE),]

    #change day on pH 7 to D7
    sub_all$day[grep("D1_7",sub_all$ind)] <- c("D7")
    # need to make replicates match. pairings are arbitrary at this point
    sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
            paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)

    Data.xtabs = xtabs(count ~ allele + day + replicate,
                   data=sub_all)
    test <- mantelhaen.test(Data.xtabs)
    d1_selection_pval[i] <- test$p.value

    #if (i%%5000 == 0){print(i)}
}

d1_selection_pval[which(is.na(d1_selection_pval))] <- 1 # bc some are invariant

print("D1-D1  cmh done")

#############
# cmh for pH selection after 7 days, comparing d7 low and control ph

d7_selection_pval <-c()

for(i in 1:nrow(mydata)){

    sub_ad <- mydata[i,grep("_DP1", colnames(mydata))]
    sub_ad <- stack(sub_ad)
    sub_ad$allele <- rep("ac1", nrow(sub_ad))
    sub_dp <- mydata[i,grep("_DP2", colnames(mydata))]
    sub_dp <- stack(sub_dp)
    sub_dp$allele <- rep("ac2", nrow(sub_dp))
    colnames(sub_dp) <- c("count", "ind", "allele")
    colnames(sub_ad) <- c("count", "ind", "allele")
    sub_all <- rbind(sub_ad, sub_dp)
    sub_all$day <- substr(sub_all$ind, 1,2)
    sub_all$replicate <- substr(sub_all$ind, 6,7)
    sub_all$pH <- substr(sub_all$ind, 4,4)

    # remove unwanted replicates
    sub_all <- sub_all[grep("D1_", sub_all$ind, invert=TRUE),]

    # need to make replicates match. pairings are arbitrary at this point
    sub_all$replicate <- ave(paste(sub_all$pH, sub_all$allele, sep=":"),
            paste(sub_all$pH, sub_all$allele, sep=":"), FUN=seq_along)

    Data.xtabs = xtabs(count ~ allele + pH + replicate,
                   data=sub_all)
    test <- mantelhaen.test(Data.xtabs)
    d7_selection_pval[i] <- test$p.value

    #if (i%%5000 == 0){print(i)}
    #ftable(Data.xtabs)
}
d7_selection_pval[which(is.na(d7_selection_pval))] <- 1 # bc some are invariant

print("D7-D7 cmh done")

########
# cmh for control selection, only include control ph indivs

control_selection_pval <-c()

for(i in 1:nrow(mydata)){

    sub_ad <- mydata[i,grep("_DP1", colnames(mydata))]
    sub_ad <- stack(sub_ad)
    sub_ad$allele <- rep("ac1", nrow(sub_ad))
    sub_dp <- mydata[i,grep("_DP2", colnames(mydata))]
    sub_dp <- stack(sub_dp)
    sub_dp$allele <- rep("ac2", nrow(sub_dp))
    colnames(sub_dp) <- c("count", "ind", "allele")
    colnames(sub_ad) <- c("count", "ind", "allele")
    sub_all <- rbind(sub_ad, sub_dp)
    sub_all$day <- substr(sub_all$ind, 1,2)
    sub_all$replicate <- substr(sub_all$ind, 6,7)
    sub_all <- sub_all[grep("D7_7_22_AD|D7_7_22_DP", sub_all$ind, invert=TRUE),] # so groups are balanced
    sub_all$pH <- substr(sub_all$ind, 4,4)
    sub_all <- sub_all[grep("_8_", sub_all$ind),]
    # need to make replicates match. pairings are arbitrary
    sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"),
            paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
    sub_all <- sub_all[grep("4", sub_all$replicate, invert=TRUE),] # subset to 3 reps to match pH treatment

    Data.xtabs = xtabs(count ~ allele + day + replicate,
                   data=sub_all)
    test <- mantelhaen.test(Data.xtabs)
    control_selection_pval[i] <- test$p.value

    #if (i%%5000 == 0){print(i)}
    #ftable(Data.xtabs)
}

control_selection_pval[which(is.na(control_selection_pval))] <- 1 # bc some are invariant

print("control cmh done")

######## append p values these to mydata ############

control_selection_pval <- control_selection_pval
mydata$pH_selection_pval <- pH_selection_pval
mydata$d1_selection_pval <- d1_selection_pval
mydata$d7_selection_pval <- d7_selection_pval

# convert to q value

mydata$control_selection_qval <- qvalue(control_selection_pval)$qvalues
mydata$pH_selection_qval <- qvalue(pH_selection_pval)$qvalues
mydata$d1_selection_qval <- qvalue(d1_selection_pval)$qvalues
mydata$d7_selection_qval <- qvalue(d7_selection_pval)$qvalues

cut_off <- 0.01

mydata$pH_sig <- FALSE
mydata$pH_sig[which(mydata$pH_selection_qval < cut_off)] <- TRUE

################ save output ########################

write.table(file="~/urchin_af/analysis/cmh.out.txt", mydata, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


##############################
##
## determining sig cutoff
##
##############################

# want an fdr of 0.01. 1%
# use control p values. these are all false positives.

length(which(mydata$control_selection_qval < cut_off))
length(which(mydata$pH_selection_qval < cut_off))
length(which(mydata$d7_selection_qval < cut_off))
length(which(mydata$d1_selection_qval < cut_off))

selected <- mydata[(which(mydata$pH_selection_qval < cut_off)),]

# write bedfile of cmh sig genes

write.table(file="~/urchin_af/analysis/cmh.all.bed",
        cbind(as.character(mydata$CHROM), mydata$POS-1, mydata$POS),
        col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file="~/urchin_af/analysis/cmh.selected.bed",
        cbind(as.character(selected$CHROM), selected$POS-1, selected$POS),
        col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")

nonselected <- mydata[(which(mydata$pH_selection_pval >= cut_off)),]

write.table(file="~/urchin_af/analysis/cmh.neutral.bed",
        cbind(as.character(nonselected$CHROM), nonselected$POS-1, nonselected$POS),
        col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")
