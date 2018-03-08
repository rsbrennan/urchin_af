snp.info <- read.table("~/urchin_af/analysis/af.info.txt", header=TRUE)
snp <- read.table("~/urchin_af/analysis/af.out.txt", header=TRUE)

colnames(snp) <- gsub("OASV2_DNA_", "", colnames(snp))
colnames(snp) <- gsub("S_", "", colnames(snp))
colnames(snp) <- gsub("5_", "", colnames(snp))
colnames(snp) <- gsub("0_", "", colnames(snp))

out.split <- apply(snp, 2, function(x) strsplit(x, ":"))

dat <- data.frame(row.names=seq(from=1, to=nrow(snp),by= 1))

for(i in 1:length(names(out.split))){
    ct <- matrix(unlist(out.split[[i]]), ncol=4, byrow=TRUE)
    #
    dat[,paste(names(out.split)[i], "DPtotal", sep="_")] <- sapply(strsplit(ct[,3], ","), "[", 1)
    dat[,paste(names(out.split)[i], "DP1", sep="_")] <- sapply(strsplit(ct[,4], ","), "[", 1)
    dat[,paste(names(out.split)[i], "DP2", sep="_")] <- sapply(strsplit(ct[,4], ","), "[", 2)
}

dep <-  dat[,grep("_DPtotal", colnames(dat))]   
dep <-  dep[,grep("D1_7_07_DPtotal", colnames(dep), invert=TRUE)]   

#convert all to numeric
dep <- data.frame(sapply(dep, function(x) as.numeric(as.character(x))))


length(which(rowSums(as.data.frame(lapply(dep,function(x)x > 50 ))) ==15))


#head((apply(dep,1, function(x) x>50)))
avgdepth <- apply(dep,1,mean)
length(which(avgdepth > 50 & avgdepth < 200))

sd(avgdepth[which(avgdepth > 50)])

length(intersect(which(avgdepth > 50 & avgdepth < 264), which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15)))

#rowSums(as.data.frame(lapply(dep,function(x)x>50))) >=14

png("~/urchin_af/figures/depth.sub.redo.hist.png", res=300, height=7, width=7, units="in")
hist(apply(dep,1,mean), breaks=500, col="grey")

dev.off()

hist(log10(avgdepth[which(avgdepth > 50)]), breaks=500, col="grey")

##
## look at effect of depth on maf. to try to get at cutoffs for depth
##

# af for looking at effect of depth on af

out <- data.frame(matrix(ncol=3, nrow=0))
colnames(out) <- c("rep", "depth", "af")

pop <- unique(substr(colnames(dep),1,7))

for(i in 1:length(pop)){
    # pull out DP1
    dp1 <- as.numeric(dep[,grep(paste(pop[i], "_DP1", sep=""), colnames(dep))])
    # total depth
    dptot <- as.numeric(dep[,grep(paste(pop[i], "_DPtotal", sep=""), colnames(dep))])
    #calc af
    af <- dp1/dptot
    #rep
    rep <- rep(pop[i], length(af))

    tmp <- data.frame(matrix(ncol=3, nrow=length(af)))
    colnames(tmp) <- c("rep", "depth", "af")

    #add to df
    tmp[,1] <- rep
    tmp[,2] <- dptot
    tmp[,3] <- af

    out <- rbind(out, tmp)

    print(i)

}



#Cut the maf values into bins by intervals of 0.05

low.out <- out[which(out$depth <= 200),]
hi.out <- out[which(out$depth > 200),]

low.out$bin <- cut(low.out$depth, breaks=seq(from=0,to=max(low.out$depth), by=20))
hi.out$bin <- cut(hi.out$depth, breaks=seq(from=200,to=max(hi.out$depth), by=200))

out.fq <- rbind(low.out, hi.out)
out.fq$bin <- as.character(out.fq$bin)
out.fq$bin[which(out.fq$depth >= 800)] <- c(">800")

out.fq$bin <- as.factor(out.fq$bin)

out.lev <-  c("(0,20]", "(20,40]","(40,60]",  "(60,80]", "(80,100]", "(100,120]", "(120,140]",
             "(140,160]", "(160,180]","(180,200]","(200,400]", "(400,600]","(600,800]",">800")
out.fq$bin<- factor(out.fq$bin, levels = out.lev)

out.fq <- out.fq[which(out.fq$depth != 0),]

# transform to maf

out.fq$af_minor <- (sapply(out.fq$af,function(x)  
          ifelse(x > 0.5, (1-x), x)))

colnames(out.fq) <- c("rep", "depth", "af", "bin", "af_minor")


library(ggplot2)

png("~/urchin_af/figures/depth.af.png", res=300, height=6, width=10, units="in")
ggplot(out.fq, aes((bin), af_minor))+
    geom_boxplot(fill="grey")+
    theme_bw()

dev.off()


#########
##
## depth calc
##
#########

#########
###

dep <-  dat[,grep("_DPtotal", colnames(dat))]   
dep <-  dep[,grep("D1_7_07_DPtotal", colnames(dep), invert=TRUE)]   

#convert all to numeric
dep <- data.frame(sapply(dep, function(x) as.numeric(as.character(x))))

#head((apply(dep,1, function(x) x>50)))
avgdepth <- apply(dep,1,mean)
#length(which(rowSums(as.data.frame(lapply(dep,function(x)x>10))) >=15))
length(which(avgdepth > 50 & avgdepth < 200))
#rowSums(as.data.frame(lapply(dep,function(x)x>50))) >=14

png("~/urchin_af/figures/depth.hist.png", res=300, height=6, width=10, units="in")
hist(avgdepth, breaks=500, col="grey")

dev.off()


png("~/urchin_af/figures/depth.sub.hist.png", res=300, height=6, width=10, units="in")
hist(avgdepth, breaks=500, col="grey", xlim=c(0,200))

dev.off()

#upper limit cutoff

mean(avgdepth[which(avgdepth > 50)])




###############
###
### output positions to include
###
###############

keep <- snp.info[intersect(which(avgdepth > 50 & avgdepth < 264), which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15)),]

write.table(file="~/urchin_af/variants/keep.snps.txt", as.data.frame(cbind(as.character(keep$CHROM), keep$POS)), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")




#############################
##
## calc allele freq to calc shifts. save output for other analyses
##
#############################


#save positions to keep

keep <- intersect(which(avgdepth > 50 & avgdepth < 264), which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15))


#### calc AF for each rep

dep <-  dat[,grep("D1_7_07", colnames(dat), invert=TRUE)]   

dep.keep <- dep[keep,]

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

dep <-  dat[,grep("D1_7_07", colnames(dat), invert=TRUE)]   

dep.keep <- dep[keep,]

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

keep.info <- snp.info[keep,1:2]


write.table(file="~/urchin_af/data/allele.freq_all.txt",cbind(keep.info, dep.keep, af, af.mean), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
