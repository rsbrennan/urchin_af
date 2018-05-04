# Fig_04_afchangecomp.R

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

tiff("~/urchin_af/figures/Fig_06_afchangecomp.tiff", height=100, width=200, units="mm", res=300)
par(mfrow = c(1, 2), mar=c(3, 3, 1.7, 1), mgp=c(3, 1, 0), las=0)
plot(0,type='n', xlim=c(0,.39), ylim=c(0,.39),
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

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.6, tcl=-0.2)
title(ylab=expression(paste(Delta," allele frequency pH 7.5")), line=1.5, cex.lab=0.7)
title(xlab=expression(paste(Delta," allele frequency pH 8.0")), line=1.5, cex.lab=0.7)

legend("topleft", c("pH 7.5 selected",
                    "pH 8.0 selected",
                    "Overlapping selected"),
    horiz = FALSE, inset = c(0, 0),
    col = c("firebrick3", "royalblue3","darkorchid4"),
    pt.bg = c("firebrick3", "royalblue3","darkorchid4"),
    pt.cex=1, cex=0.7, pch=c(21,22,24))


mtext(text="A",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)

######
##
## maf plot
##
######
 #sampling from entire dist

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", stringsAsFactors=FALSE, header=TRUE)

nrep <- 1000
nsamp <- length(snp.sel_75)

bs <- matrix(nrow=nsamp, ncol=nrep)
avg_rep <- c()
for(i in 1:nrep){
    bs[,i] <- sample(mydata$af_out, size=nsamp, replace=FALSE)
    avg_rep <- c(avg_rep, bs[,i])
}

out <- list()
for (i in 1: ncol(bs)){

    out[[i]] <- data.frame(x=density(bs[,i], bw = 0.05)$x, 
        y= density(bs[,i], bw = 0.05)$y, 
        bin = cut(density(bs[,i], bw = 0.05)$x, breaks=seq(from=-0.3, to=1.3, by=0.01)))
}

out.new <- do.call(rbind, out)

densities.qtiles <- out.new %>% 
group_by(bin) %>%
  summarise(q05 = quantile(y, 0.025),
            q50 = quantile(y, 0.5),
            q95 = quantile(y, 0.975))

densities.qtiles$x <- seq(from=-0.155, to=1.135, by=0.01)


avg_perm <- c()
for (i in 1: ncol(bs)){
    avg_perm <- c(avg_perm, bs[,i])
}



plot(density(0:1), ylim=c(0,4),xlim=c(0,1), lwd=0,
    main="",
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.6, tcl=-0.2)
title(xlab="Starting allele frequency", line=1.5, cex.lab=0.7)

polygon(x=c(densities.qtiles$x,rev(densities.qtiles$x)),
    y=c(densities.qtiles$q05,rev(densities.qtiles$q95)),
    col=alpha("black", alpha=0.2),border=NA)

lines(densities.qtiles$x,densities.qtiles$q50, col="black", lwd=3 )
#lines(density(mydata$af_out), col=alpha("black", 1), lwd=3)
lines(density(snp.sel_75, bw=0.05), col=alpha("firebrick3", 1), lwd=3)
lines(density(snp.sel_80,bw=0.05), col=alpha("royalblue3", 1), lwd=3)
lines(density(snp.sel_both,bw=0.05), col=alpha("darkorchid2", 1), lwd=3)

legend("topright", c("permuted", "pH 7.5: selected", "pH 8.0: selected", "Overlapping selected"),
    col=c("black", "firebrick3", "royalblue3", "darkorchid2"), lty=1,
    cex=0.7, lwd=2)


mtext(text="B",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.14*diff(par("usr")[1:2]), outer=FALSE)


dev.off()


