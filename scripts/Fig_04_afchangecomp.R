library(scales)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
cut_off <- 0.001

selected_7 <- mydata[which(mydata$pH7_selection_qval < cut_off & mydata$pH8_selection_qval >= cut_off),]
selected_8 <- mydata[which(mydata$pH8_selection_qval < cut_off & mydata$pH7_selection_qval >= cut_off),]
selected_both <- mydata[which(mydata$pH8_selection_qval < cut_off & mydata$pH7_selection_qval < cut_off),]

selected_7_D7_8 <- selected_7[,grep("D7_8_mean", colnames(selected_7))]
selected_7_D7_7 <- selected_7[,grep("D7_7_mean", colnames(selected_7))]
selected_7_D1_8 <- selected_7[,grep("D1_8_mean", colnames(selected_7))]

selected_8_D7_8 <- selected_8[,grep("D7_8_mean", colnames(selected_8))]
selected_8_D7_7 <- selected_8[,grep("D7_7_mean", colnames(selected_8))]
selected_8_D1_8 <- selected_8[,grep("D1_8_mean", colnames(selected_8))]

selected_both_D7_8 <- selected_both[,grep("D7_8_mean", colnames(selected_both))]
selected_both_D7_7 <- selected_both[,grep("D7_7_mean", colnames(selected_both))]
selected_both_D1_8 <- selected_both[,grep("D1_8_mean", colnames(selected_both))]

selected_all <- unique(rbind(selected_7,selected_8,selected_both))
selected_all_D7_8 <- selected_all[,grep("D7_8_mean", colnames(selected_all))]
selected_all_D7_7 <- selected_all[,grep("D7_7_mean", colnames(selected_all))]
selected_all_D1_8 <- selected_all[,grep("D1_8_mean", colnames(selected_all))]


d7_8 <- abs(selected_all_D1_8 - selected_all_D7_8)
d7_7 <- abs(selected_all_D1_8 - selected_all_D7_7)

d7_8_both <- abs(selected_both_D1_8 - selected_both_D7_8)
d7_7_both <- abs(selected_both_D1_8 - selected_both_D7_7)

d7_8_s7 <- abs(selected_7_D1_8  - selected_7_D7_8)
d7_7_s7 <- abs(selected_7_D1_8  - selected_7_D7_7)

d7_8_s8 <- abs(selected_8_D1_8  - selected_8_D7_8)
d7_7_s8 <- abs(selected_8_D1_8  - selected_8_D7_7)


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

png("~/urchin_af/figures/Fig_06_afchangecomp.png", height=100, width=120, units="mm", res=300)

par(mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
par(fig = c(0,1,0,1)) # this sets location of first plot
plot(0,type='n', xlim=c(0,.39), ylim=c(0,.39),
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")
box(which="plot")
points(x=d7_8_s7, y=d7_7_s7, pch=21, col=alpha("firebrick3", 0.2),
    bg = alpha("firebrick3", 0.2), cex=0.8)
points(x=d7_8_s8, y=d7_7_s8, pch=21, col=alpha("royalblue3", 0.2),
    bg = alpha("royalblue3", 0.2), cex=0.8)
points(x=d7_8_both, y=d7_7_both, pch=21, col=alpha("darkorchid4", 0.7),
    bg = alpha("darkorchid4", 0.7), cex=0.8)

abline(0, 1, col="black", lty=2, lwd=2.2)

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.6, tcl=-0.2)
title(ylab=expression(paste(Delta," allele frequency pH 7.5")), line=1.5, cex.lab=0.7)
title(xlab=expression(paste(Delta," allele frequency pH 8.0")), line=1.5, cex.lab=0.7)

legend("topleft", c("pH 7.5 significant",
                    "pH 8.0 significant",
                    "Overlapping selected"),
    horiz = FALSE, inset = c(0, 0), pch = c(19, 19, 19),
    col = c("firebrick3", "royalblue3","darkorchid4"), pt.cex=1, cex=0.7)

dev.off()
