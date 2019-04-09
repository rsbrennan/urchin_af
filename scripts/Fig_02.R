# Fig_02

# makes fig 2 a and b.

# pca and scatter of af changes


library(pcadapt)
library(dplyr)
library(scales)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)

dat <- mydata[,grep("_af", colnames(mydata))]
dat_t <- t(dat)

filename <- read.pcadapt(dat_t,type="pool")
x <- pcadapt(filename,K=5) # calc actual pca


poplist.names <- row.names(dat_t)
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


# get proportion of total variance explained:
x$singular.values[1]/sum(x$singular.values)*100
x$singular.values[2]/sum(x$singular.values)*100


# black circle Day 1 8.0
# red circle Day 1 7.5
# black triangle Day 7 8.0
# red triangle Day 7 7.5

sp <- c(21,23,22)
bg.col <-  c("firebrick3", "gray48","royalblue3")

tiff("~/urchin_af/figures/Fig_02.tiff", height=85, width=170, units="mm", res=300)
par(mfrow = c(1, 2), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)

plot(y=pc_sum$PC2, x=pc_sum$PC1,
    pch=sp[as.factor(pc_sum$pH_day)],
    cex=1.6, col="black",
    bg=bg.col[as.factor(pc_sum$pH_day)],
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    #ylim=c(-0.085, 0.11),
    #xlim=c(-0.086, 0.125),
    xaxt="n",yaxt="n"
)

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)

title(xlab="PC1: 22.0%",  line=1.5, cex.lab=0.9)
title(ylab="PC2: 20.7%",  line=1.5, cex.lab=0.9)

legend("bottomright", legend=c(expression('T'[0]), "pH 7.5", "pH 8.0" ),
    pt.bg=c("gray31", "firebrick3", "royalblue3") ,
    pt.cex=1.3, cex=0.6, pch=c(23,21,22),
    bg="white")

mtext(text=bquote(paste('(',italic('a'),')')),
              side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)

# plot the af change results

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)

cut_off <- 0.05/(9828)

mydata$D1_8_maf <- (sapply(mydata$D1_8_af,function(x)
          ifelse(x > 0.5, (1-x), x)))

selected_7 <- mydata[which(mydata$pH7_selection_pval < cut_off & mydata$pH8_selection_pval >= cut_off),]
selected_8 <- mydata[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval >= cut_off),]
selected_both <- mydata[which(mydata$pH8_selection_pval < cut_off & mydata$pH7_selection_pval < cut_off),]
neutral <- mydata[which(mydata$pH8_selection_pval > cut_off & mydata$pH7_selection_pval > cut_off),]

# pull out allele freqs
selected_7_D7_8 <- selected_7$D7_8_af
selected_7_D7_7 <- selected_7$D7_7_af
selected_7_D1_8 <- selected_7$D1_8_af

selected_8_D7_8 <- selected_8$D7_8_af
selected_8_D7_7 <- selected_8$D7_7_af
selected_8_D1_8 <- selected_8$D1_8_af

selected_both_D7_8 <- selected_both$D7_8_af
selected_both_D7_7 <- selected_both$D7_7_af
selected_both_D1_8 <- selected_both$D1_8_af

selected_all <- unique(rbind(selected_7,selected_8,selected_both))
selected_all_D7_8 <- selected_all$D7_8_af
selected_all_D7_7 <- selected_all$D7_7_af
selected_all_D1_8 <- selected_all$D1_8_af

d7_8 <- (selected_all_D7_8 - selected_all_D1_8)
d7_7 <- (selected_all_D7_7 - selected_all_D1_8)

d7_8_both <- (selected_both_D7_8 - selected_both_D1_8)
d7_7_both <- (selected_both_D7_7 - selected_both_D1_8)

d7_8_s7 <- (selected_7_D7_8 - selected_7_D1_8)
d7_7_s7 <- (selected_7_D7_7 - selected_7_D1_8)

d7_8_s8 <- (selected_8_D7_8 - selected_8_D1_8)
d7_7_s8 <- (selected_8_D7_7 - selected_8_D1_8)


cor.test(d7_7, d7_8,
         method = "pearson",
         conf.level = 0.95)
cor.test(d7_7_s7, d7_8_s7,
         method = "pearson",
         conf.level = 0.95)
cor.test(d7_7_s8, d7_8_s8,
         method = "pearson",
         conf.level = 0.95)

cor.test(d7_7_both, d7_8_both,
         method = "pearson",
         conf.level = 0.95)


# corr genome wide:

d7_8_all <- (neutral$D7_8_af - neutral$D1_8_af)
d7_7_all <- (neutral$D7_7_af - neutral$D1_8_af)

cor.test(d7_8_all, d7_7_all)




#########
#
# plot figure 4, folded allele freq change
#
#########

#tiff("~/urchin_af/figures/Fig_04_folded_afchangecomp.tiff", height=85, width=170, units="mm", res=300)
plot(0,type='n', xlim=c(-0.03,.39), ylim=c(-0.03,.39),
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")
box(which="plot")
points(x=d7_8_s7, y=d7_7_s7, col=alpha("firebrick3", 0.2),
    bg = alpha("firebrick3", 0.2), cex=0.8, pch=21)
points(x=d7_8_s8, y=d7_7_s8, col=alpha("royalblue3", 0.2),
    bg = alpha("royalblue3", 0.2),  cex=0.8, pch=22)
points(x=d7_8_both, y=d7_7_both, col=alpha("darkorchid4", 0.7),
    bg = alpha("darkorchid4", 0.7), cex=0.8, pch=24)

abline(0, 1, col="black", lty=2, lwd=2.2)
abline(h=0, lty=2, col="grey")
abline(v=0, lty=2, col="grey")
axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)
title(ylab=expression(paste(Delta," allele frequency pH 7.5")), line=1.5, cex.lab=0.9)
title(xlab=expression(paste(Delta," allele frequency pH 8.0")), line=1.5, cex.lab=0.9)

legend("topleft", c("pH 7.5 selected",
                    "pH 8.0 selected",
                    "Overlapping selected"),
    horiz = FALSE, inset = c(0, 0),
    pt.bg = c("firebrick3", "royalblue3","darkorchid4"),
    pt.cex=1.3, cex=0.6, pch=c(21,22,24),
    bg="white")

mtext(text=bquote(paste('(',italic('b'),')')),
              side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)

dev.off()


