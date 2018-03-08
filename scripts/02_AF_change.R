library(dplyr)
library(reshape)
library(ggplot2)
library(qvalue)
library(scales)
library(qqman)


# obtain allele freqs:

command1 <- "zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^##' | cut -f 1-8 | sed 's/#CHROM/CHROM/g'"
snp.info <- read.table(pipe(command1), header=TRUE)
command2 <- "zcat ~/urchin_af/variants/urchin_final.vcf.gz | grep -v '^##' | cut -f 10-"
snp <- read.table(pipe(command2), header=TRUE)

colnames(snp) <- gsub("OASV2_DNA_", "", colnames(snp))
colnames(snp) <- gsub("S_", "", colnames(snp))
colnames(snp) <- gsub("5_", "", colnames(snp))
colnames(snp) <- gsub("0_", "", colnames(snp))

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
dat <- data.frame(sapply(dat, function(x) as.numeric(as.character(x))))

#### calc AF for each rep
dep.keep <- dat

# make empty data frame to hold af
af <- data.frame(matrix(ncol=15, nrow=nrow(dep.keep)))
pop <- unique(substr(colnames(dep),1,7))
colnames(af) <- pop

for(i in 1:ncol(af)){
    # pull out DP1
    dp1 <- as.numeric(dep.keep[,grep(paste(pop[i], "_DP1", sep=""), colnames(dep.keep))])
    # total depth
    dptot <- as.numeric(dep.keep[,grep(paste(pop[i], "_DPtotal", sep=""), colnames(dep.keep))])
    #calc af
    a.freq <- dp1/dptot
    # add to af
    af[,grep(paste(pop[i]), colnames(af))] <- a.freq
    print(i)
}

colnames(af) <- paste(pop, "_af", sep="")

# calculate mean of each group
#### calc AF for each rep

# calculate mean of each group
gp <- c("D1_7", "D1_8", "D7_7", "D7_8")
af.mean <-  data.frame(matrix(ncol=length(gp), nrow=nrow(dep.keep)))
colnames(af.mean) <- gp

for (i in 1:length(gp)){
    sub.gp <- dep.keep[,grep(gp[i], colnames(dep.keep))]
    #subset pops
    pop <- unique(substr(colnames(sub.gp),1,7))
    #make tmp matrix
    tmp.out <-  data.frame(matrix(ncol=length(pop), nrow=nrow(dep.keep)))
    colnames(tmp.out) <- pop
    for (subpop in 1:length(pop)){
        # pull out DP1
        dp1 <- as.numeric(sub.gp[,grep(paste(pop[subpop], "_DP1", sep=""), colnames(sub.gp))])
        # total depth
        dptot <- as.numeric(sub.gp[,grep(paste(pop[subpop], "_DPtotal", sep=""), colnames(sub.gp))])
        #calc af
        a.freq <- dp1/dptot
        tmp.out[,pop[subpop]] <- a.freq
    }
    af.mean[gp[i]] <- apply(tmp.out,1,mean)
}
colnames(af.mean) <- c("D1_7_mean", "D1_8_mean", "D7_7_mean", "D7_8_mean")

# add snp info, etc to this.
keep.info <- snp.info[,1:2]

write.table(file="~/urchin_af/data/allele.freq.txt",cbind(keep.info, dep.keep, af, af.mean), 
  row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


######
# actual af change analysis


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

mydata$control_selection_pval <- control_selection_pval
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
