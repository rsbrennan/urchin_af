
# subsampling down to 3 samples per rep.
# removing D7_8_17 because pca outlier
# want to see if enrichment in p values for low pH treatment still holds

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
#	count <- sub_dp$depth - sub_ad$count
#	sub_ac2 <- data.frame(count=count, ind=sub_dp$ind, allele=sub_dp$allele) 
	sub_all <- rbind(sub_ad, sub_dp)
	sub_all$day <- substr(sub_all$ind, 1,2)
	sub_all$replicate <- substr(sub_all$ind, 6,7)
	sub_all$pH <- substr(sub_all$ind, 4,4)

	# remove unwanted replicates
	sub_all <- sub_all[grep("D1_7_", sub_all$ind, invert=TRUE),]
	sub_all <- sub_all[grep("D7_8_", sub_all$ind, invert=TRUE),]
	# randomly remove one replicate from each
	rm_d1_8 <- sample(unique(substr(sub_all$ind[grep("D1_8",sub_all$ind)], 1,7)), 1)
	sub_all <- sub_all[grep(rm_d1_8, sub_all$ind, invert=TRUE),]
	rm_d7_7 <- sample(unique(substr(sub_all$ind[grep("D7_7",sub_all$ind)], 1,7)), 1)
	sub_all <- sub_all[grep(rm_d7_7, sub_all$ind, invert=TRUE),]
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

	# randomly remove one replicate from each to get down to 3 samps. also remove D7_8_17
	rm_d1_8 <- sample(unique(substr(sub_all$ind[grep("D1_8",sub_all$ind)], 1,7)), 1)
	sub_all <- sub_all[grep(rm_d1_8, sub_all$ind, invert=TRUE),]
	sub_all <- sub_all[grep("D7_8_17", sub_all$ind, invert=TRUE),]

	# need to make replicates match. pairings are arbitrary
	sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"), 
			paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)
	sub_all <- sub_all[grep("4", sub_all$replicate, invert=TRUE),] # subset to 3 reps to match pH treatment

	Data.xtabs = xtabs(count ~ allele + day + replicate, 
                   data=sub_all)
	test <- mantelhaen.test(Data.xtabs)
	control_selection_pval[i] <- test$p.value

	if (i%%5000 == 0){print(i)}
	#ftable(Data.xtabs)
}

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
	# randomly remove one replicate from each to get down to 3 samps. also remove D7_8_17
	rm_d1_8 <- sample(unique(substr(sub_all$ind[grep("D1_8",sub_all$ind)], 1,7)), 1)
	sub_all <- sub_all[grep(rm_d1_8, sub_all$ind, invert=TRUE),]

	#change day on pH 7 to D7
	sub_all$day[grep("D1_7",sub_all$ind)] <- c("D7")
	# need to make replicates match. pairings are arbitrary at this point
	sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"), 
			paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)

	Data.xtabs = xtabs(count ~ allele + day + replicate, 
                   data=sub_all)
	test <- mantelhaen.test(Data.xtabs)
	d1_selection_pval[i] <- test$p.value

	if (i%%5000 == 0){print(i)}
}

######## append p values these to mydata ############
mydata$control_selection_pval <- control_selection_pval
mydata$pH_selection_pval <- pH_selection_pval
mydata$d1_selection_pval <-d1_selection_pval

################ save output ########################
write.table(file="~/urchin_af/analysis/cmh_3samp.out.txt", mydata, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

mydata <- read.table("~/urchin_af/analysis/cmh_3samp.out.txt", header=TRUE)

png("~/urchin_af/figures/pval_hist_3samps.png", res=300, height=7, width=10, units="in")
par(mfrow = c(1, 1))

plot(density(mydata$pH_selection_pval, na.rm=TRUE),col="red",
	main="CMH p-value- 3 samples",lwd=2.5)
lines(density(mydata$control_selection_pval, na.rm=TRUE), col="black", lwd=2.5)
legend("topright", c("D1-D7 pH 7.5", "D1-D7 pH 8.0"), pch=19, col=c("red", "black", "blue"))

dev.off()

