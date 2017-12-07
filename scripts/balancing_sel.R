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
