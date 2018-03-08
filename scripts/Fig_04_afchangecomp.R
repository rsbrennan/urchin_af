library(scales)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
cut_off <- 0.01

selected <- mydata[(which(mydata$pH_selection_qval < cut_off)),]
selected <- selected[(which(selected$control_selection_qval >= cut_off)),]

D1_8 <- selected[,grep("D1_8_mean", colnames(selected))]

D1_7 <- selected[,grep("D1_7_mean", colnames(selected))]

D7_7 <- selected[,grep("D7_7_mean", colnames(selected))]
D7_8 <- selected[,grep("D7_8_mean", colnames(selected))]

d1_d7 <- abs(D1_8-D7_7)
d1_d1 <- abs(D1_8-D1_7)
d1_d7_ctr <- abs(D1_8-D7_8)

cor.test(d1_d7, d1_d1, 
         method = "pearson",
         conf.level = 0.95)

cor.test(d1_d7, d1_d7_ctr, 
         method = "pearson",
         conf.level = 0.95)

cor.test(d1_d7, d1_d1, 
         method = "spearman",
         conf.level = 0.95)

cor.test(d1_d7, d1_d1, 
         method = "kendall",
         conf.level = 0.95)


reg1 <- lm(d1_d7~d1_d1)

png("~/urchin_af/figures/Fig_06_afchangecomp.tiff", height=100, width=100, units="mm", res=300)

par(mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
par(fig = c(0,1,0,1)) # this sets location of first plot
plot(0,type='n', xlim=c(0,0.35), ylim=c(0,0.35),
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")  
box(which="plot")
points(x=d1_d1, y=d1_d7, pch=19, col = alpha("black", 0.5), cex=0.7)
#abline(reg1, col="firebrick3", lwd=2)
abline(0, 1, col="firebrick3", lty=2, lwd=2)

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.6, tcl=-0.2)
title(ylab=expression(paste(Delta," allele frequency: ", "D1 pH 8.0 vs. D7 pH 7.5")), line=1.5, cex.lab=0.7)
title(xlab=expression(paste( Delta," allele frequency: ", "D1 pH 8.0 vs. D1 pH 7.5")), line=1.5, cex.lab=0.7)

#legend("topleft", c("1:1", "Observed"), xpd = TRUE, 
#   horiz = FALSE, inset = c(0, 0), 
#   lty = c(2, 1), col = c("firebrick3", "firebrick3"), lwd=1.8, cex=0.6)

# add rectangle around inset fig
rect(xleft=0.17, xright=0.365, ybottom=-0.1, ytop=0.16)

# add inset histogram
par(fig = c(0.46,.98, 0.065, 0.58), new = T)  # location of 2nd plot. fist 2 nums are x locaiton, 2nd are y loc

hist(d1_d1,  col = alpha("black", 0.5), freq=FALSE,
    ylab="", xlab="",
    xaxt="n",yaxt="n",
    main="",
    xlim=c(0,.3))   
hist(d1_d7,  col = alpha("firebrick3", 0.5), freq=FALSE, add=T)
axis(1, mgp=c(0.5, - 0.25, 0), cex.axis=0.5, tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(0.5, 0.09, 0), cex.axis=0.5, tcl=-0.2)

title(ylab="Density", cex.lab=0.7, line=0.55, cex.lab=0.5)
title(xlab="Avg. change in allele frequency", line=0.3, cex.lab=0.5)
legend("topright", c("D1 pH 8.0 vs. D1 pH 7.5", "D1 pH 8.0 vs. D7 pH 7.5"), 
    horiz = FALSE, inset = c(0, 0), 
    bty = "n", pch = c(15, 15), col = c("black", "firebrick3"), pt.cex=0.7, cex=0.5)


# add legend. basically overlaying new empty plot
#par(fig = c(0.13, 1, 0, 1), oma = c(0, 0, .15, 0), mar = c(0, 0, 0, 0), new = TRUE)
#plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("topleft", c("1:1", "Observed"), xpd = TRUE, 
#   horiz = FALSE, inset = c(0, 0), 
#   bty = "n", lty = c(2, 1), col = c("red", "red"), lwd=1.8, cex=0.6)

#text(x=0.5, y=.95, paste("y = ",round(coef(reg1)[2], 2), "x", " + ", round(coef(reg1)[1], 2), sep=""), cex=0.7, xpd = TRUE)
#text(x=-0.1, y=.95, paste("R-squared:", round(summary(reg1)$adj.r.squared, 2)), cex=0.7)

dev.off()



########
# separate plots, for presentation fig
########


pdf("~/urchin_af/figures/Fig_06_1.pdf", height=5, width=5)

par(mfrow=c(1,1))

par(mar=c(5, 5, 1.7, 1), mgp=c(3, 1, 0), las=0)
plot(0,type='n', xlim=c(0,0.3), ylim=c(0,0.3),
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")  
box(which="plot")
points(x=d1_d1, y=d1_d7, pch=19, col = alpha("black", 0.5), cex=1)
#abline(reg1, col="firebrick3", lwd=2)
abline(0, 1, col="firebrick3", lty=2, lwd=2)

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.8,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.8, tcl=-0.2)
title(ylab=expression(paste(Delta," allele frequency: ", "D1 pH 8.0 vs. D7 pH 7.5")), line=2.3, cex.lab=1)
title(xlab=expression(paste( Delta," allele frequency: ", "D1 pH 8.0 vs. D1 pH 7.5")), line=2.3, cex.lab=1)

#legend("bottomright", c("1:1"), xpd = TRUE, 
#   horiz = FALSE, inset = c(0, 0), 
#   lty = c(2, 1), col = c("firebrick3", "firebrick3"), lwd=1.8, cex=1)

dev.off()
# add rectangle around inset fig
#rect(xleft=0.17, xright=0.365, ybottom=-0.1, ytop=0.16)

# add inset histogram
pdf("~/urchin_af/figures/Fig_06_2.pdf", height=5, width=5)

hist(d1_d1,  col = alpha("black", 0.5), freq=FALSE,
    ylab="", xlab="",
    xaxt="n",yaxt="n",
    main="",
    xlim=c(0,.3))   
hist(d1_d7,  col = alpha("firebrick3", 0.5), freq=FALSE, add=T)
axis(1, mgp=c(0.5, 0.1, 0), cex.axis=0.8, tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(0.5, 0.1, 0), cex.axis=0.8, tcl=-0.2)

title(ylab="Density", cex.lab=1, line=2)
title(xlab="Avg. change in allele frequency", line=2, cex.lab=1)
legend("topright", c("D1 pH 8.0 vs. D1 pH 7.5", "D1 pH 8.0 vs. D7 pH 7.5"), 
    horiz = FALSE, inset = c(0, 0), 
    bty = "n", pch = c(15, 15), col = c("black", "firebrick3"), pt.cex=1, cex=1)

dev.off()

