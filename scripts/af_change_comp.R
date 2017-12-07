
# generate plot of AF change 
# compare D1-D7 comparison vs. D1-D1 comparison


mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
cut_off <- quantile(mydata$control_selection_pval, 0.01, na.rm=TRUE)

selected <- mydata[(which(mydata$pH_selection_pval < cut_off)),]

delta_d7 <- abs(selected$D1_8_mean-selected$D7_7_mean)
delta_d1 <- abs(selected$D1_8_mean-selected$D1_7_mean)

png("~/urchin_af/figures/af_change_comp.png", res=300, height=7, width=7, units="in")

par(mar=c(6, 5.7, 4.1, 2.1))
plot(delta_d7, delta_d1, pch=19,
    xlim=c(0,.31), ylim=c(0,0.31),
    ylab="Mean change in allele frequency\n day1 pH 8.0 vs. day 1 pH 7.5",
    xlab="",
    cex.lab=1.2)

title(xlab="Mean change in allele frequency\n day1 pH 8.0 vs. day 7 pH 7.5", 
    line=3.8, cex.lab=1.2,)

abline(0, 1, lwd=3, lty=2, col='red') 
abline(lm(delta_d1~delta_d7), col="red", lwd=3)

legend("topleft", c("1:1 relationship", "observed relationship"), lty=2:1, col=c("red", "red"), cex=1.1, lwd=3)

dev.off()
