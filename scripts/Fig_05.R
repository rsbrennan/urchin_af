# Fig_05.R

library(qqman)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
mydata$CHR <- as.numeric(gsub("Scaffold", "", mydata$CHROM))
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")


cut_off <- 0.001

#allele frequencies to plot panel
#calc stderror of mean af estimates

gps <- c("D1_7", "D7_7", "D7_8")

af_dat <- mydata[,grep("_af", colnames(mydata))]

mydata$D7_7_delta <- abs(mydata$D1_8_mean-mydata$D7_7_mean)
mydata$D7_8_delta <- abs(mydata$D1_8_mean-mydata$D7_8_mean)

new.dat <- subset(mydata, CHR > 541 & CHR < 543)
new.dat <- subset(new.dat, POS > 34350 & POS < 34510)

d=data.frame(CHR=new.dat$CHR, BP=new.dat$POS,
            P=new.dat$pH7_selection_qval, P_8= new.dat$pH8_selection_qval,
            SNP=new.dat$SNP,
            D7_7_delta = new.dat$D7_7_delta, D7_8_delta=new.dat$D7_8_delta)

    d <- d[order(d$CHR, d$BP), ]
    d$logp <- d$P
    d$pos=NA
    d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency
    nchr = length(unique(d$CHR))
    if (nchr==1) { ## For a single chromosome
        d$pos=d$BP
        xlabel = paste('Chromosome',unique(d$CHR),'position')
    } else { ## For multiple chromosomes
        lastbase=0
        ticks=NULL
        for (i in unique(d$index)) {
            if (i==1) {
                d[d$index==i, ]$pos=d[d$index==i, ]$BP
            } else {
        lastbase = lastbase +max(d[d$index==(i-1),"BP"])
        d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
        d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase
            }
        }
    #ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
    labs <- unique(d$CHR)
    }

    # Initialize plot
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)

# define plot parameters:
    def_args <- list(xaxt='n', yaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                     xlim=c(xmin,xmax), ylim=c(0,16),
                     xlab="", ylab="")

tiff("~/urchin_af/figures/Fig_05.tiff", res=300, height=109, width=85, units="mm")
#dev.new(height=5, width=3.9)

par(mfrow = c(2, 1), mar=c(1, 4, 3, 1), mgp=c(3, 1, 0), las=0)

plot(0,type='n', xlim=c(34350,34510), ylim=c(0,.2),
    main="",
    ylab="",
    xlab="",
    cex.lab=1, cex.axis=0.7,
    xaxt="n",yaxt="n")

axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2, labels=FALSE) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)

lines((d$pos), d$D7_7_delta, col="firebrick3")
points(d$pos, d$D7_7_delta, col="black",
    bg = "firebrick3", cex=1, pch=21)
lines((d$pos), d$D7_8_delta, col="royalblue3")
points(d$pos, d$D7_8_delta, col="black",
    bg = "royalblue3", cex=1, pch=22)

mtext(side=2,expression(paste(Delta, " allele frequency")), line=2.2, cex=0.8)
mtext(side=2, expression("vs. T"[0]), line=1.4, cex=0.8)

# add label
mtext(text="A",
        side=3, line=0.3,
             cex=1.5,
            at=par("usr")[1]-0.2*diff(par("usr")[1:2]))

legend("top",c("pH 7.5", "pH 8.0"), pch = c(21,22),
    col="black", pt.bg=c("firebrick3", "royalblue3"),
    horiz = TRUE, inset = c(0, -0.2), xpd=TRUE, cex=0.9,
    bty = "n", pt.cex=1.4)



###########
#pval plot:
###########
par(mar=c(3, 4, 1, 1), mgp=c(3, 1, 0), las=0)

plot(0,type='n', xlim=c(34350,34510), ylim=c(0,4.3),
    main="",
    ylab="",
    xlab="",
    cex.lab=1, cex.axis=0.7,
    xaxt="n",yaxt="n")
axis(1, mgp=c(1.8, .2, 0), cex.axis=0.7,tcl=-0.2) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7, tcl=-0.2)

abline(h=-log10(cut_off), lty=2, col="black", lwd=1.5)

lines((d$pos), -log10(d$P), col="firebrick3")
points(d$pos, -log10(d$P), col="black",
    bg = "firebrick3", cex=1, pch=21)
lines((d$pos), -log10(d$P_8), col="royalblue3")
points(d$pos, -log10(d$P_8), col="black",
    bg = "royalblue3", cex=1, pch=22)

title(ylab=expression(-log[10](italic(q))), line=2, cex.lab=0.8)

title(xlab="Scaffold542: position in BP", line=1.9, cex.lab=0.8)
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7)


# add label
mtext(text="B",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.2*diff(par("usr")[1:2]))

dev.off()
