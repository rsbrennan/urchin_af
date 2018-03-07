


snp.info <- read.table("~/urchin_af/analysis/af.info.txt", header=TRUE)
snp <- read.table("~/urchin_af/analysis/af.out.txt", header=TRUE)

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

dep <-  dat[,grep("_DPtotal", colnames(dat))]
# remove bad sample w/ low coverage
dep <-  dep[,grep("D1_7_07_DPtotal", colnames(dep), invert=TRUE)]

#convert all to numeric
dep <- data.frame(sapply(dep, function(x) as.numeric(as.character(x))))

print("number of site with > 50x coverage in all samples")
length(which(rowSums(as.data.frame(lapply(dep,function(x)x > 50 ))) ==15))

avgdepth <- apply(dep,1,mean)
print("number of sites with avg depth > 50, < 264")
length(which(avgdepth > 50 & avgdepth < 264))

print("sd of sites with avg depth > 50")
sd(avgdepth[which(avgdepth > 50)])

print("mean depth of sites with avg depth > 50")
mean(avgdepth[which(avgdepth > 50 )])

print("filtering based on sites w/ mean depth> 50, max depth < mean depth * 3, or 88*3, 264")

length(intersect(which(avgdepth > 50 & avgdepth < 264), which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15)))

#########
##
## depth calc
##
#########

dep <-  dat[,grep("_DPtotal", colnames(dat))]
dep <-  dep[,grep("D1_7_07_DPtotal", colnames(dep), invert=TRUE)]

#convert all to numeric
dep <- data.frame(sapply(dep, function(x) as.numeric(as.character(x))))

avgdepth <- apply(dep,1,mean)

png("~/urchin_af/figures/depth.hist.png", res=300, height=6, width=10, units="in")
hist(avgdepth, breaks=500, col="grey")

dev.off()

png("~/urchin_af/figures/depth.sub.hist.png", res=300, height=6, width=10, units="in")
hist(avgdepth, breaks=500, col="grey", xlim=c(0,200))
dev.off()

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

