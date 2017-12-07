# cmh plots


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
## plot all responsive loci
##
##############################

selected <- mydata[(which(mydata$pH_selection_pval < cut_off)),]

nonselected <- mydata[(which(mydata$pH_selection_pval >= cut_off)),]


# thin to those that are showing response in treatment
selected <- selected[(which(mydata$pH_selection_pval < cut_off)),]
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
## manhattan plots
##
##############################



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

##############################
##
## glm plots
##
##############################

dat_glm <- read.table("~/urchin_af/analysis/glm.out.txt", header=TRUE)

dat <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

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

png("~/urchin_af/figures/glm_interaction_hist.png", res=300, height=7, width=7, units="in")
hist(dat_glm$glm_pval, col='grey', breaks=40, main="glm day*pH interaction",
	xlab="day*pH p-value")
dev.off()

png("~/urchin_af/figures/glm_ph_hist.png", res=300, height=7, width=7, units="in")
hist(dat_glm$glm_ph_pval, col='grey', breaks=40, main="glm pH",
	xlab="pH p-value")
dev.off()




