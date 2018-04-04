library(pcadapt)
library(dplyr)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

dat <- mydata[,grep("_af", colnames(mydata))]
dat_t <- t(dat)

filename <- read.pcadapt(dat_t,type="pool",local.env = TRUE,pop.sizes = rep(50, nrow(dat_t)))
x <- pcadapt(filename,K=5) # calc actual pca

# pop names
poplist.int <- c(rep(1,50),rep(2,50),rep(3,50))
# With names
poplist.names <- c(rep(row.names(dat_t)[1],50),
                rep(row.names(dat_t)[2],50),
                rep(row.names(dat_t)[3],50),
                rep(row.names(dat_t)[4],50),
                rep(row.names(dat_t)[5],50),
                rep(row.names(dat_t)[6],50),
                rep(row.names(dat_t)[7],50),
                rep(row.names(dat_t)[8],50),
                rep(row.names(dat_t)[9],50),
                rep(row.names(dat_t)[10],50),
                rep(row.names(dat_t)[11],50),
                rep(row.names(dat_t)[12],50))

pop <- substr(poplist.names, 1,4)
ph <- substr(poplist.names, 4,4)
day <- substr(poplist.names, 1,2)
num <- substr(poplist.names, 1,7)
pH_day <- paste(ph, day, sep="_")
dat.p <- data.frame(pop=pop , ph = ph, day= day, num=num,pH_day =pH_day, PC1 = x$scores[,1],  PC2= x$scores[,2])

pc1 <-group_by(dat.p, num) %>% summarize(PC1 = mean(PC1))
pc2 <-group_by(dat.p, num) %>% summarize(PC2 = mean(PC2))

pc_sum <- as.data.frame(cbind(pc1, pc2))
pc_sum$pH_day <- paste(substr(as.character(pc_sum$num), 4,4),substr(as.character(pc_sum$num), 1,2), sep="_")

# jitter D7 ph8 point so we can see
pc_sum$PC1[10] <- (pc_sum$PC1[10]- pc_sum$PC1[10])

# flip each pc
pc_sum$PC1 <- pc_sum$PC1*(-1)
pc_sum$PC2 <- pc_sum$PC2*(-1)

# get proportion of total variance explained:
x$singular.values[1]/sum(x$singular.values)*100
x$singular.values[2]/sum(x$singular.values)*100


# black circle Day 1 8.0
# red circle Day 1 7.5
# black triangle Day 7 8.0
# red triangle Day 7 7.5


sp <- c(24,21,22)
bg.col <-  c("firebrick3", "gray48","royalblue3")

tiff("~/urchin_af/figures/Fig_03_pca.tiff", res=300, height=85, width=85, units="mm")

#dev.new(width=3.35, height=3.35, units="mm")
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(y=pc_sum$PC2, x=pc_sum$PC1,
    pch=sp[as.factor(pc_sum$pH_day)],
    cex=1.6, col="black",
    bg=bg.col[as.factor(pc_sum$pH_day)],
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    ylim=c(-0.085, 0.11),
    xlim=c(-0.086, 0.125),
    xaxt="n",yaxt="n"
)

axis(1, mgp=c(2, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(2, .5, 0), cex.axis=0.7)

title(xlab="PC1: 22.0%", line=2, cex.lab=1.05)
title(ylab="PC2: 20.7%", line=2, cex.lab=1.05)

legend("topright", pch=c(21,22,24),
    pt.bg=c("gray31", "royalblue3", "firebrick3") ,
    legend=c(expression('T'[0]), "pH 8.0", "pH 7.5" ),
    pt.cex=1.6, cex=0.85)

dev.off()
