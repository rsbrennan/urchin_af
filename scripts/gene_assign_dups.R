dat <- read.table("~/urchin_af/analysis/cmh.all.genes.txt", header=FALSE)
dat$SNP <- paste(dat$V1, dat$V3, sep=":")
snp <- unique(dat$SNP)
out <- as.data.frame(matrix(ncol=ncol(dat), nrow=0))
z <- 0
for (i in 1:length(snp)) {
    a <- dat[which(dat$SNP == snp[i]),]
    if (length(unique(a$V14)) > 1){
        z <- z+1
        print(z)
        genes <- unique(a$V14)
        for (g in 1:length(genes)){
            g_tmp <- which(a$V14 == genes[g])
            out <- rbind(out, a[g_tmp[1],])
        }
    }
    else{
        out <- rbind(out, a[1,])
    }
    if (i%%1000 == 0) {print(i)}
}

write.table(file="~/urchin_af/analysis/cmh.all_dups_in.genes.txt", out, col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
