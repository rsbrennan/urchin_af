


# read in data
snp.info <- read.table("~/urchin_af/analysis/af.info.txt", header=TRUE)
snp <- read.table("~/urchin_af/analysis/af.out.txt", header=TRUE)

# make names a bit less annoying
colnames(snp) <- gsub("OASV2_DNA_", "", colnames(snp))
colnames(snp) <- gsub("S_", "", colnames(snp))
colnames(snp) <- gsub("5_", "", colnames(snp))
colnames(snp) <- gsub("0_", "", colnames(snp))

#split af
out.split <- apply(snp, 2, function(x) strsplit(x, ":"))

#make new data frame to save output
dat <- data.frame(row.names=seq(from=1, to=nrow(snp),by= 1))

# get depths and fill dat data frame
for(i in 1:length(names(out.split))){
    ct <- matrix(unlist(out.split[[i]]), ncol=4, byrow=TRUE)
    #
    dat[,paste(names(out.split)[i], "DPtotal", sep="_")] <- sapply(strsplit(ct[,3], ","), "[", 1)
    dat[,paste(names(out.split)[i], "DP1", sep="_")] <- sapply(strsplit(ct[,4], ","), "[", 1)
    dat[,paste(names(out.split)[i], "DP2", sep="_")] <- sapply(strsplit(ct[,4], ","), "[", 2)
}

dep <-  dat[,grep("_DPtotal", colnames(dat))]
#convert all to numeric
dep <- data.frame(sapply(dep, function(x) as.numeric(as.character(x))))

avgdepth <- apply(dep,1,mean)

mean(avgdepth)
