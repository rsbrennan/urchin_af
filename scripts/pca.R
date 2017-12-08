library(pcadapt)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

dat <- mydata[,grep("_af", colnames(mydata))] 
dat_t <- t(dat)

filename <- read.pcadapt(dat_t,type="pool",local.env = TRUE,pop.sizes = rep(50, nrow(dat_t)))
x <- pcadapt(filename,K=20) # calc actual pca

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
                rep(row.names(dat_t)[12],50),
                rep(row.names(dat_t)[13],50),
                rep(row.names(dat_t)[14],50),
                rep(row.names(dat_t)[15],50))   

pop <- substr(poplist.names, 1,4)
ph <- substr(poplist.names, 4,4)
day <- substr(poplist.names, 1,2)
num <- substr(poplist.names, 1,7)
pH_day <- paste(ph, day, sep="_")
dat.p <- data.frame(pop=pop , ph = ph, day= day, num=num,pH_day =pH_day, PC1 = x$scores[,1],  PC2= x$scores[,2])

p <- ggplot(dat.p , aes(x= PC1, y=PC2)) + 
    geom_point(aes(fill = as.factor(dat.p$num), shape = as.factor(pH_day)),size = 4, color="black") + 
    scale_shape_manual(values=c(21,22, 24,25)) + 
    theme_bw()+
    xlab("PC1: 7.2%") + #x$singular.values[1]/sum(x$singular.values) 
    ylab("PC2: 6.6%") +
    guides(fill=FALSE)
    #theme(text = element_text(size=18))
#   labs(color='pH', shape='pH')+
 #  ggtitle(paste("locus", loci, sep="-"))
png("~/urchin_af/figures/pca.png", res=300, height=7, width=7, units="in")
p
dev.off()

#par(mfrow = c(1, 1))
plot(x,option="scores",pop=poplist.names)
plot(x,option="screeplot")

# rerun pca without outlier sample: D7_8_17

dat_t <- t(dat)
dat_r <- dat_t[grep("D7_8_17", rownames(dat_t), invert=TRUE),]

filename <- read.pcadapt(dat_r,type="pool",local.env = TRUE,pop.sizes = rep(50, nrow(dat_r)))
x <- pcadapt(filename,K=20) # calc actual pca

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
                rep(row.names(dat_t)[12],50),
                rep(row.names(dat_t)[13],50),
                rep(row.names(dat_t)[14],50))   

pop <- substr(poplist.names, 1,4)
ph <- substr(poplist.names, 4,4)
day <- substr(poplist.names, 1,2)
num <- substr(poplist.names, 1,7)
pH_day <- paste(ph, day, sep="_")
dat.p <- data.frame(pop=pop , ph = ph, day= day, num=num,pH_day =pH_day, PC1 = x$scores[,1],  PC2= x$scores[,2])

p <- ggplot(dat.p , aes(x= PC1, y=PC2)) + 
    geom_point(aes(fill = as.factor(dat.p$num), shape = as.factor(pH_day)),size = 4, color="black") + 
    scale_shape_manual(values=c(21,22, 24,25)) + 
    theme_bw()+
    xlab("PC1: 7.1%") + #x$singular.values[1]/sum(x$singular.values) 
    ylab("PC2: 6.5%") +
    guides(fill=FALSE)
    #theme(text = element_text(size=18))
#   labs(color='pH', shape='pH')+
 #  ggtitle(paste("locus", loci, sep="-"))
png("~/urchin_af/figures/pca_outlier_rm.png", res=300, height=7, width=7, units="in")
p
dev.off()
