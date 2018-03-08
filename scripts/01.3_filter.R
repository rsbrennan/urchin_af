
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

print("number of site with > 40x coverage in all samples")
length(which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15))

avgdepth <- apply(dep,1,mean)
print("number of sites with avg depth > 50, < 264")
length(which(avgdepth >= 50 & avgdepth < 264))

print("mean depth of sites with avg depth > 50")
mean(avgdepth[which(avgdepth > 50 )])

print("filtering based on sites w/ mean depth> 50, max depth < mean depth * 3, 
    or 88*3, 264, and each pool sequenced > 40x")

length(intersect(which(avgdepth > 50 & avgdepth < 264), which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15)))

###############
###
### output positions to include
###
###############

keep <- snp.info[intersect(which(avgdepth > 50 & avgdepth < 264), which(rowSums(as.data.frame(lapply(dep,function(x)x > 40 ))) ==15)),]

write.table(file="~/urchin_af/variants/keep.snps.txt", as.data.frame(cbind(as.character(keep$CHROM), keep$POS)), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
