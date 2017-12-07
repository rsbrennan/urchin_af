
library(dplyr)
library(reshape) 
library(ggplot2)
library(qvalue)
library(scales)
library(qqman)

# analysis of sea urchin AF change. 1 generation. 

mydata <- read.table("~/urchin_af/data/allele.freq.txt", header=TRUE) 
head(mydata)

# the colnames are awful. get rid of some of the mess so they make sense

## subset data to only include af measures

dat <- mydata[,grep("_af", colnames(mydata))] 

#########
##
## identify variants only responsive to pH treatment
##
#########

# from Whole-Genome Resequencing of Experimental Populations Reveals Polygenic Basis of Egg-Size Variation in Drosophila melanogaster
	# we first identified variants that have shifted in the same directions relative to the starting 
	# population in all three replicates in each treatment and then calculated the average allele 
	# frequency change relative to the starting population across each treatment, that is, 
	# (abs(average(p1-starting population, p2-starting population, p3-starting population)), 
	# where p = LEP, SEP, and CP). To calculate allele frequency difference between the LEP and SEP, 
	# we first identified variants that have shifted in the same direction relative to the starting population 
	# in all three replicate populations of the LEP and SEP. We then calculated difference in average allele 
	# frequencies between the two treatments.

# make empty data frame
af.dir <- as.data.frame(matrix(ncol=3, nrow=nrow(dat)))
colnames(af.dir) <- c("pH_response", "lab_response", "d1_response")

for(i in 1:nrow(dat)){
	sub_7.trt <- dat[i,grep("D7_7_", colnames(dat))]
	sub_8.trt <- dat[i,grep("D7_8_", colnames(dat))]
	sub_7.ctr <- mean(t(dat[i,grep("D1_8_", colnames(dat))]))
	sub_8.ctr <- mean(t(dat[i,grep("D1_8_", colnames(dat))]))
	if(sum(sub_7.trt > sub_8.ctr) == 4 | sum(sub_7.trt < sub_8.ctr) == 4) {
		af.dir$pH_response[i] <- TRUE
	} else{
		af.dir$pH_response[i] <- FALSE
	} 
	if(sum(sub_8.trt > sub_8.ctr) == 4 | sum(sub_8.trt < sub_8.ctr) == 4) {
		af.dir$lab_response[i] <- TRUE
	} else{
		af.dir$lab_response[i] <- FALSE
	} 
	if(sum(sub_7.ctr > sub_8.ctr) == 3 | sum(sub_7.ctr < sub_8.ctr) == 3) {
		af.dir$d1_response[i] <- TRUE
	} else{
		af.dir$d1_response[i] <- FALSE
	} 
   if (i%%5000 == 0){print(i)} # printing progress
}

for(i in 1:nrow(dat)){
	sub_7.trt <- dat[i,grep("D1_7_", colnames(dat))]
	sub_8.trt <- dat[i,grep("D7_8_", colnames(dat))]
	sub_7.ctr <- mean(t(dat[i,grep("D1_7_", colnames(dat))]))
	sub_8.ctr <- mean(t(dat[i,grep("D1_8_", colnames(dat))]))
	if(sum(sub_7.trt > sub_8.ctr) == 3 | sum(sub_7.ctr < sub_8.ctr) == 3) {
		af.dir$d1_response[i] <- TRUE
	} else{
		af.dir$d1_response[i] <- FALSE
	} 
   if (i%%5000 == 0){print(i)}
}

# determine if shifts are in same direction or not
af.dir$pH_selected <-ifelse((af.dir$pH_response == TRUE & af.dir$lab_response == FALSE), TRUE, FALSE)
af.dir$lab_selected <-ifelse((af.dir$pH_response == FALSE & af.dir$lab_response == TRUE), TRUE, FALSE)
af.dir$d1_selected <-ifelse((af.dir$d1_response == TRUE & af.dir$lab_response == FALSE), TRUE, FALSE)

mydata <- cbind(mydata, af.dir)

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
#	count <- sub_dp$depth - sub_ad$count
#	sub_ac2 <- data.frame(count=count, ind=sub_dp$ind, allele=sub_dp$allele) 
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

	if (i%%5000 == 0){print(i)}
	#ftable(Data.xtabs)
}

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

######## append p values these to mydata ############

mydata$control_selection_pval <- control_selection_pval
mydata$pH_selection_pval <- pH_selection_pval
mydata$d1_selection_pval <- d1_selection_pval
mydata$d7_selection_pval <- d7_selection_pval

cut_off <- quantile(mydata$control_selection_pval, 0.01, na.rm=TRUE)

mydata$pH_sig <- FALSE
mydata$pH_sig[which(mydata$pH_selection_pval < cut_off)] <- TRUE

################ save output ########################

write.table(file="~/urchin_af/analysis/cmh.out.txt", mydata, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

##### Plotting ######
mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

png("~/urchin_af/figures/pval_hist.png", res=300, height=7, width=10, units="in")

par(mfrow = c(1, 1))

plot(density(mydata$pH_selection_pval, na.rm=TRUE),col="red",
	main="CMH p-value density plots",lwd=2.5)
lines(density(mydata$control_selection_pval, na.rm=TRUE), col="black", lwd=2.5)
lines(density(mydata$d7_selection_pval, na.rm=TRUE), col="blue", lwd=2.5)
legend("topright", c("D1-D7 pH 7.5", "D1-D7 pH 8.0", "D7_7.5-D7_8.0"), pch=19, col=c("red", "black", "blue"))

dev.off()

##### plot individual loci ########

### read in data:
mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

test <- mydata[which(mydata$pH_selection_pval < 0.001),]
length(which(test$lab_response ==TRUE))
length(which(test$pH_response ==TRUE))
#length(which(test$d1_response ==TRUE))

test <- test[(which(test$pH_response ==TRUE)),]

test.af <- test[,grep("_af", colnames(test))] 

loci <- grep("Scaffold542:34445", test$SNP)

#11, 19 is good

x <- stack(test.af[loci,]);
xx <- cbind(x$values, colsplit(x$ind, "_", names = c("x", "y", "z", "zz")));
names(xx) <- c("allelefreq", "day", "pH", "replicate", "ig");
#xx$day <- as.numeric(gsub("D", "", xx$day))
p <- ggplot(xx, aes(x= day, y=allelefreq, shape = as.factor(pH), color = as.factor(pH), group=as.factor(pH))) + 
	stat_smooth(se=FALSE,method="lm") +
	stat_summary(fun.data="mean_se", mapping=aes(group=list(interaction(day, as.factor(pH)))), 
		size=1.5,show.legend=FALSE, color="black", position = position_dodge(width = .5))+
	geom_point(position=position_jitter(w=0.15,h=0), size = 4) + 
	theme_bw() +
	theme(text = element_text(size=18))+
 	labs(color='pH', shape='pH')+
 	ggtitle(paste("locus", loci, sep="-"))
p
ggsave("~/urchin_af/figures/indiv.task2_2.png", plot=p)


##############################
##
## determining sig cutoff
##
##############################

# want an fdr of 0.01. 1%
# use control p values. these are all false positives. 

cut_off <- quantile(mydata$control_selection_pval, 0.01, na.rm=TRUE)
length(which(mydata$control_selection_pval < cut_off))
length(which(mydata$pH_selection_pval < cut_off))
length(which(mydata$d7_selection_pval < cut_off))


##############################
##
## plot all responsive loci
##
##############################

selected <- mydata[(which(mydata$pH_selection_pval < cut_off)),]

# write bedfile of cmh sig genes

write.table(file="~/urchin_af/analysis/cmh.selected.bed", 
            cbind(as.character(selected$CHROM), selected$POS-1, selected$POS), 
            col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")

nonselected <- mydata[(which(mydata$pH_selection_pval >= cut_off)),]

write.table(file="~/urchin_af/analysis/cmh.neutral.bed", 
            cbind(as.character(nonselected$CHROM), nonselected$POS-1, nonselected$POS), 
            col.names=FALSE, row.names=FALSE, quote=FALSE,sep="\t")

# thin to those that are showing response in treatment
selected <- selected[(which(selected$pH_response == TRUE)),]
length(which(mydata$pH_selection_pval < cut_off))
length(which(mydata$control_selection_pval < cut_off))
length(which(mydata$d1_selection_pval < cut_off))
selected.af <- selected[,grep("_mean", colnames(selected))] 

selected.t <- selected.af[,grep("D1_8_mean|D7_7_mean", colnames(selected.af))] 
selected.treat <- data.frame(D7_7_mean=selected.t$D7_7_mean, D1_8_mean=selected.t$D1_8_mean)
mean(abs(selected.treat[,1]-selected.treat[,2]))

max(abs(selected.treat[,1]-selected.treat[,2]))

mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")
highlight <- mydata$SNP[which(mydata$pH_response == TRUE & mydata$pH_selection_pval < cut_off)]

png("~/urchin_af/figures/pH_response_ext.png", res=300, height=7, width=10, units="in")
par(mfrow = c(1, 2))

plot(x=c(0,0), y=c(1,1),
	ylim=c(0,1), xlim=c(8.25,6.75), type="n",
	xlab="pH", ylab="allele freq", xaxt = "n")
	axis(1, at=c(7,8), labels=c("treatment", "control"),
		main="Sig. AF shifts in pH only: down")

for (i in 1:nrow(selected.treat)){
	if(selected.treat[i,1]< selected.treat[i,2]){
		lines(x=c(7, 8), y=selected.treat[i,],
			col = alpha("blue", 0.2))
	}
}

# add avg response
avg.resp <- mean(abs(selected.treat[,1]-selected.treat[,2]))

lines(x=c(7, 8), y=c(mean(selected.treat$D1_8_mean)-avg.resp, mean(selected.treat$D1_8_mean)),
		col = alpha("black", 1), lwd=4, lty=2)
lines(x=c(7, 8), y=c(mean(selected.treat$D1_8_mean)-max(abs(selected.treat[,1]-selected.treat[,2])), mean(selected.treat$D1_8_mean)),
		col = alpha("black", 1), lwd=4, lty=1)

plot(x=c(0,0), y=c(1,1),
	ylim=c(0,1), xlim=c(8.25,6.75), type="n",
	xlab="pH", ylab="allele freq", xaxt = "n")
	axis(1, at=c(7,8), labels=c("treatment", "control"))

for (i in 1:nrow(selected.treat)){
	if(selected.treat[i,1]> selected.treat[i,2]){
		lines(x=c(7, 8), y=selected.treat[i,],
			col = alpha("red", 0.2))
	} 
}

lines(x=c(7, 8), y=c(mean(selected.treat$D1_8_mean)+avg.resp,mean(selected.treat$D1_8_mean)),
		col = alpha("black", 1), lwd=4, lty=2)

#extreme
lines(x=c(7, 8), 
	y=c((mean(selected.treat$D1_8_mean)-0.3)+max(abs(selected.treat[,1]-selected.treat[,2])),(mean(selected.treat$D1_8_mean)-0.3)),
		col = alpha("black", 1), lwd=4, lty=1)

dev.off()


##############################
##
## plot only pH response
##
##############################

selected <- mydata[(which(mydata$pH_selection_pval < cut_off)),]

# thin to those that are showing response in treatment
selected <- selected[(which(selected$pH_selected == TRUE)),]

selected.af <- selected[,grep("_mean", colnames(selected))] 

selected.t <- selected.af[,grep("D7_7_|D1_8_", colnames(selected.af))] 
selected.treat<- data.frame(D7_7_mean=selected.t$D7_7_mean, D1_8_mean=selected.t$D1_8_mean)
png("~/urchin_af/figures/selected.png", res=300, height=7, width=10, units="in")
par(mfrow = c(1, 2))

plot(x=c(0,0), y=c(1,1),
	ylim=c(0,1), xlim=c(8.25,6.75), type="n",
	xlab="pH", ylab="allele freq", xaxt = "n")
	axis(1, at=c(7,8), labels=c("treatment", "control"),
		main="Sig. AF shifts in pH only: down")

for (i in 1:nrow(selected.treat)){
	if(selected.treat[i,1]< selected.treat[i,2]){
		lines(x=c(7, 8), y=selected.treat[i,],
			col = alpha("blue", 0.2))
	}
}

plot(x=c(0,0), y=c(1,1),
	ylim=c(0,1), xlim=c(8.25,6.75), type="n",
	xlab="pH", ylab="allele freq", xaxt = "n")
	axis(1, at=c(7,8), labels=c("treatment", "control"),
		main="Sig. AF shifts in pH only: up")

for (i in 1:nrow(selected.treat)){
	if(selected.treat[i,1]> selected.treat[i,2]){
		lines(x=c(7, 8), y=selected.treat[i,],
			col = alpha("red", 0.2))
	} 
}
dev.off()


### calc mean






### manhattan plots

mydata$CHR <- as.numeric(gsub("Scaffold", "", mydata$CHROM))
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

png("~/urchin_af/figures/genome_wide_manhattan.png", res=300, height=5, width=10, units="in")

highlight <- mydata$SNP[which(mydata$pH_response == TRUE & mydata$pH_selection_pval < cut_off)]
manhattan(subset(mydata), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold",
	main="")

dev.off()

#highlight acid potassium channel

png("~/urchin_af/figures/genome_wide_manhattan-1.png", res=300, height=5, width=10, units="in")

highlight <- c("Scaffold542:34370", "Scaffold542:34445","Scaffold542:34450","Scaffold542:34467","Scaffold542:34469" )

manhattan(subset(mydata), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold",
	main="")

dev.off()

png("~/urchin_af/figures/manhattan_cand.png", res=300, height=5, width=10, units="in")
highlight <- c("Scaffold542:34370", "Scaffold542:34445","Scaffold542:34450","Scaffold542:34467","Scaffold542:34469" )
manhattan(subset(mydata, CHR > 538 & CHR < 543), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold")

dev.off()


highlight <- mydata$SNP[which(mydata$pH_response == TRUE & mydata$pH_selection_pval < cut_off)]

png("~/urchin_af/figures/manhattan_1.png", res=300, height=5, width=10, units="in")
manhattan(subset(mydata, CHR > 701 & CHR < 715), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold")

dev.off()

png("~/urchin_af/figures/manhattan_2.png", res=300, height=5, width=10, units="in")

manhattan(subset(mydata, CHR > 65 & CHR < 68), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold")

dev.off()

png("~/urchin_af/figures/manhattan_3.png", res=300, height=5, width=10, units="in")
manhattan(subset(mydata, CHR > 12 & CHR < 40), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold")

dev.off()

png("~/urchin_af/figures/manhattan_4.png", res=300, height=5, width=10, units="in")
manhattan(subset(mydata, CHR > 663 & CHR < 672), chr="CHR", bp="POS", p="pH_selection_pval",
	genomewideline = -log10(cut_off),
	suggestiveline = FALSE, 
	highlight= highlight,
	ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold")

dev.off()










########################################################################################
############################################
##
## glm
## 
##
############################################
########################################################################################

mydata <- read.table("~/urchin_af/data/allele.freq.txt", header=TRUE) 

glm_pval <-c()
glm_ph_pval <- c()

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

	#sub_all <- sub_all[grep("D1_7_", sub_all$ind, invert=TRUE),]

	sub_all$day <- substr(sub_all$ind, 1,2)
	sub_all$replicate <- substr(sub_all$ind, 6,7)
	sub_all$pH <- substr(sub_all$ind, 4,4)

	sub_all$replicate <- ave(paste(sub_all$day, sub_all$allele, sep=":"), 
			paste(sub_all$day, sub_all$allele, sep=":"), FUN=seq_along)

	sub_all$mer <- substr(sub_all$ind, 1, 7)
	d1 <- subset(sub_all, allele== "ac1")
	d2 <- subset(sub_all, allele== "ac2")
	dat_in <- merge(d1, d2, by="mer")

	res <- glm(cbind(dat_in$count.x, dat_in$count.y) ~ pH.x+day.x+pH.x:day.x,
		family="quasibinomial", data=dat_in)

	# save results
	n_rows<-nrow(summary(res)$coefficients)
	glm_ph_pval[i] <- summary(res)$coefficients[2,4]
	glm_pval[i] <- summary(res)$coefficients[n_rows,4]

	if (i%%5000 == 0){print(i)}
	#ftable(Data.xtabs)
}
 
mydata$glm_pval <- glm_pval
mydata$glm_ph_pval <- glm_ph_pval

qobj <- qvalue(p = dat_glm$glm_ph_pval)
mydata$glm_qval <- qobj$qvalues

mydata$CHR <- as.numeric(gsub("Scaffold", "", mydata$CHROM))
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

png("~/urchin_af/figures/glm_pval_manhattan.png", res=300, height=5, width=10, units="in")

#highlight <- mydata$SNP[which(mydata$pH_response == TRUE & mydata$pH_selection_pval < cut_off)]
manhattan(mydata, chr="CHR", bp="POS", p="glm_pval",
	genomewideline = -log10(0.05),
	suggestiveline = FALSE, 
	#highlight= highlight,
	ylim=c(0, (-log10(min(mydata$glm_pval, na.rm=TRUE)))*1.1),
	xlab="Scaffold",
	main="",
	logp=TRUE)

dev.off()

write.table(file="~/urchin_af/analysis/glm.out.txt", mydata, col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

dat_glm <- read.table("~/urchin_af/analysis/glm.out.txt", header=TRUE)

dat <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

png("~/urchin_af/figures/glm_interaction_hist.png", res=300, height=7, width=7, units="in")
hist(dat_glm$glm_pval, col='grey', breaks=40, main="glm day*pH interaction",
	xlab="day*pH p-value")
dev.off()

png("~/urchin_af/figures/glm_ph_hist.png", res=300, height=7, width=7, units="in")
hist(dat_glm$glm_ph_pval, col='grey', breaks=40, main="glm pH",
	xlab="pH p-value")
dev.off()


############################################
##
## look at starting af
## balancing selection
##
############################################

# compare pH selection vs control

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

qobj <- qvalue(p = mydata$control_selection_pval)
qval <- qobj$qvalues
q_p <- data.frame(pval=mydata$control_selection_pval, qval=qobj$qvalues)

length(which(q_p$qval < 0.01))
cut_off <- max(q_p$pval[which(q_p$qval < 0.01)])


selected <- mydata[(which(mydata$pH_selection_pval < cut_off)),]
# thin to those that are showing response in treatment
selected <- selected[(which(selected$pH_selected == TRUE)),]


dat <- selected[,grep("_af", colnames(selected))] 
dat <- dat[,grep("D1_8", colnames(dat))] 

sel.sfs = as.data.frame(lapply(dat,function(x)  
          ifelse(x > 0.5, (1-x), x)))

all.dat <- mydata[,grep("_af", colnames(mydata))] 
all.dat <- all.dat[,grep("D1_8", colnames(all.dat))] 

all.sfs = as.data.frame(lapply(all.dat,function(x)  
          ifelse(x > 0.5, (1-x), x)))

random <- matrix(nrow=nrow(selected), ncol=4)

perm.mean <- c()
perm.median <- c()
perm.75 <- c()
perm.90 <- c()
perm.95 <- c()
perm.10 <- c()
# permutation of groups. compare starting af summary stats of observed to these
# this is just randomly assigning replicates to groups, then pulling out the same
	# number of snps as under "selection" in pH

for (i in 1: 5000){
	random <- matrix(nrow=nrow(selected), ncol=4)
	random <- all.sfs[sample(nrow(all.sfs), size=nrow(random), replace=FALSE),]

	perm.mean[i] <- mean(apply(random,1,mean))
	perm.median[i] <- median(apply(random,1,mean))
	perm.75[i] <- quantile(apply(random,1,mean), 0.75)
	perm.90[i] <- quantile(apply(random,1,mean), 0.90)
	perm.95[i] <- quantile(apply(random,1,mean), 0.95)
	perm.10[i] <- quantile(apply(random,1,mean), 0.10)
	if (i%%500 == 0){print(i)}
}

control <- mydata[(which(mydata$control_selection_pval < cut_off)),]
control <- control[(which(control$lab_selected == TRUE)),]


dat <- control[,grep("_af", colnames(control))] 
dat <- dat[,grep("D1_8", colnames(dat))] 

ctr.sfs = as.data.frame(lapply(dat,function(x)  
          ifelse(x > 0.5, (1-x), x)))


png("~/urchin_af/figures/pH_permutations.png", res=300, height=15, width=10, units="in")
par(mfrow = c(3, 2))

hist(perm.mean, col="grey", breaks=100, xlim=c(0.17, 0.25))
abline(v=mean(apply(sel.sfs,1,mean)), col="red", lwd=3)
abline(v=mean(apply(ctr.sfs,1,mean)), col="green", lty=2, lwd=3)

hist(perm.median, col="grey", breaks=100, xlim=c(0.12, 0.22))
abline(v=median(apply(sel.sfs,1,mean)), col="red", lwd=3)
abline(v=median(apply(ctr.sfs,1,mean)), col="green", lty=2, lwd=3)

hist(perm.75, col="grey", breaks=100,
	xlim=c(0.2, .35))
abline(v=quantile(apply(sel.sfs,1,mean), 0.75), col="red", lwd=3)
abline(v=quantile(apply(ctr.sfs,1,mean), 0.75), col="green", lty=2, lwd=3)

hist(perm.90, col="grey", breaks=100,
	xlim=c(0.35, quantile(apply(sel.sfs,1,mean), 0.90)*1.3))
abline(v=quantile(apply(sel.sfs,1,mean), 0.90), col="red", lwd=3)
abline(v=quantile(apply(ctr.sfs,1,mean), 0.90), col="green", lty=2, lwd=3)

hist(perm.95, col="grey", breaks=100,
	xlim=c(0.4, quantile(apply(sel.sfs,1,mean), 0.95)*1.1))
abline(v=quantile(apply(sel.sfs,1,mean), 0.95), col="red", lwd=3)
abline(v=quantile(apply(ctr.sfs,1,mean), 0.95), col="green", lty=2, lwd=3)

hist(perm.10, col="grey", breaks=100,
	xlim=c(0, 0.09))
abline(v=quantile(apply(sel.sfs,1,mean), 0.1), col="red", lwd=3)
abline(v=quantile(apply(ctr.sfs,1,mean), 0.1), col="green", lty=2, lwd=3)

dev.off()

### plot his of control and treat maf

png("~/urchin_af/figures/maf_hist.png", res=300, height=7, width=7, units="in")
par(mfrow = c(1, 1))

hist(apply(sel.sfs,1,mean), breaks=60, freq=FALSE, col = alpha("red", 0.4), main="pH selection (red) vs. control (black)")
hist(apply(ctr.sfs,1,mean), breaks=60, freq=FALSE, col = alpha("black", 0.4), add=T)
#hist(d1_selection_pval, breaks=80, freq=FALSE, col = alpha("orange", 0.4), add=T)

dev.off()













################################################
######
## permute to pull out false positives and compare sfs
######
################################################

# permute samples. pull out "responsive" loci. calc summary stats
# try to run in parallel

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-6) #not to overload your computer
registerDoParallel(cl)

#loop
my_results_par <- foreach(perm_rep = 1:500, .combine = rbind) %dopar%
    {  
#    sink("~/monitor.txt",append=TRUE)
	cat(paste("Starting iteration",perm_rep,"\n"), 
       file="~/log.monitor.txt", append=TRUE)

	ad <- mydata[,grep("_AD", colnames(mydata))]	
	ad.dp <- mydata[,grep("_DP", colnames(mydata))]
	shuf <- sample(seq(1:ncol(ad)))
	ad <- ad[,shuf]
	ad.dp <- ad.dp[,shuf]

	ad.sub <- ad[,1:8]
	colnames(ad.sub) <- c("D1_1_ac","D1_2_ac","D1_3_ac","D1_4_ac","D7_1_ac","D7_2_ac","D7_3_ac","D7_4_ac")
	ad.dp.sub <- ad.dp[,1:8]
	colnames(ad.dp.sub) <- c("D1_1_dp","D1_2_dp","D1_3_dp","D1_4_dp","D7_1_dp","D7_2_dp","D7_3_dp","D7_4_dp")
	for(i in 1:nrow(ad.sub)){

		sub_ad <- stack(ad.sub[i,])
		sub_ad$allele <- rep("ac1", nrow(sub_ad))
		sub_dp <- stack(ad.dp.sub[i,])
		sub_dp$allele <- rep("ac2", nrow(sub_dp))
		colnames(sub_dp) <- c("depth", "ind", "allele")
		colnames(sub_ad) <- c("count", "ind", "allele")
		count <- sub_dp$depth - sub_ad$count
		sub_ac2 <- data.frame(count=count, ind=sub_dp$ind, allele=sub_dp$allele) 
		sub_all <- rbind(sub_ad, sub_ac2)
		sub_all$day <- substr(sub_all$ind, 1,2)
		sub_all$replicate <- substr(sub_all$ind, 4,4)
		Data.xtabs = xtabs(count ~ allele + day + replicate, 
    	               data=sub_all)
		test <- mantelhaen.test(Data.xtabs)
		control_selection_pval[i] <- test$p.value
	
		#if (i%%20000 == 0){cat(paste("On subiteration",perm_rep, ":", i,"\n"))}

		#ftable(Data.xtabs)
	}

	# pull out sfs of sig results

	cut_off <- quantile(control_selection_pval, 0.001, na.rm=TRUE)
	out <- cbind(ad.sub,ad.dp.sub )
	selected <- out[(which(control_selection_pval < cut_off)),]
	fre <- as.data.frame(matrix(nrow=nrow(out), ncol=4))
		fre[,1] <- out$D1_1_ac/out$D1_1_dp
		fre[,2] <- out$D1_2_ac/out$D1_2_dp
		fre[,3] <- out$D1_3_ac/out$D1_3_dp
		fre[,4] <- out$D1_4_ac/out$D1_4_dp

	perm.sfs = as.data.frame(lapply(fre,function(x)  
          ifelse(x > 0.5, (1-x), x)))

	perm.mean <- mean(apply(perm.sfs,1,mean))
	perm.median <- median(apply(perm.sfs,1,mean))
	perm.75 <- quantile(apply(perm.sfs,1,mean), 0.75)
	perm.90 <- quantile(apply(perm.sfs,1,mean), 0.90)
	perm.95 <- quantile(apply(perm.sfs,1,mean), 0.95)
	cbind(perm.mean,perm.median,perm.75,perm.90,perm.95)
}

write.table(my_results_par, file="~/urchin_af/analysis/permutation.txt",  col.names=TRUE, quote=FALSE, sep="\t")

#stop cluster
stopCluster(cl)

my_results_par <- as.data.frame(my_results_par)


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


########################################
##
## LD calculations
##
########################################


# MBE: Genomics of Parallel Experimental Evolution in Drosophila

need to use LDx



