# pi columns are: chrom, window position (mean val), number of snps in window, fraction of site with sufficient cov, pi measure for window

D1_8 <- read.table("~/urchin_af/variants/D1_8.pi", header=FALSE,stringsAsFactors=FALSE)
D1_8 <- D1_8[which(D1_8$V5 != "na"),]
colnames(D1_8) <- c("CHR", "POS", "n_snp", "per_cov", "pi")
mean(as.numeric(D1_8$pi))

# pi: 0.01582153

D7_7 <- read.table("~/urchin_af/variants/D7_7.pi", header=FALSE, stringsAsFactors=FALSE)
D7_7 <- D7_7[which(D7_7$V5 != "na"),]
colnames(D7_7) <- c("CHR", "POS", "n_snp", "per_cov", "pi")

mean(as.numeric(D7_7$pi))

D7_8 <- read.table("~/urchin_af/variants/D7_8.pi", header=FALSE, stringsAsFactors=FALSE)
D7_8 <- D7_8[which(D7_8$V5 != "na"),]
colnames(D7_8) <- c("CHR", "POS", "n_snp", "per_cov", "pi")

mean(as.numeric(D7_8$pi))

# fairly different number per group
nrow(D1_8)
nrow(D7_7)
nrow(D7_8)

D1_8$snp <- paste(D1_8$CHR, D1_8$POS, sep=":")
D7_7$snp <- paste(D7_7$CHR, D7_7$POS, sep=":")
D7_8$snp <- paste(D7_8$CHR, D7_8$POS, sep=":")

d1 <- merge(D1_8, D7_7, by="snp")
d2 <- merge(d1, D7_8, by="snp")

mean(as.numeric(d2$pi.x)) # D1_8
mean(as.numeric(d2$pi.y)) # D7_7
mean(as.numeric(d2$pi)) # D7_8

wilcox.test(as.numeric(d2$pi.x), as.numeric(d2$pi.y), correct=TRUE)
wilcox.test(as.numeric(d2$pi.x), as.numeric(d2$pi), correct=TRUE)
wilcox.test(as.numeric(d2$pi.y), as.numeric(d2$pi), correct=TRUE)
