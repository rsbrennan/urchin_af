library(qqman)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
mydata$CHR <- as.numeric(gsub("Scaffold", "", mydata$CHROM))
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

cut_off <- quantile(mydata$control_selection_pval, 0.01, na.rm=TRUE)

#allele frequencies to plot panel
#calc stderror of mean af estimates

gps <- c("D1_7", "D7_7", "D7_8")

af_dat <- mydata[,grep("_af", colnames(mydata))] 

af_se <- matrix(ncol=4, nrow=nrow(af_dat))

for (i in 1:length(gps)){
    af_tmp <- af_dat[,grep(gps[i], colnames(af_dat))] 
    sd_tmp <- apply(af_tmp, 1, sd)
    af_se[,i] <- sd_tmp/ncol(af_tmp)
}



mydata$D7_7_delta <- abs(mydata$D1_8_mean-mydata$D7_7_mean)
mydata$D1_7_delta <- abs(mydata$D1_8_mean-mydata$D1_7_mean)
mydata$D7_8_delta <- abs(mydata$D1_8_mean-mydata$D7_8_mean)

new.dat <- subset(mydata, CHR > 541 & CHR < 543)

d=data.frame(CHR=new.dat$CHR, BP=new.dat$POS, P=new.dat$pH_selection_pval, SNP=new.dat$SNP,
            D7_7_delta = new.dat$D7_7_delta,D1_7_delta = new.dat$D1_7_delta, D7_8_delta=new.dat$D7_8_delta) 
    
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

d2 <- d
#highlight <- subset(d, P < cut_off)

png("~/urchin_af/figures/Fig_07_manhattan.png",width=6.65, height= 4, units="in", res=300)

m <- cbind(c(1,2,3), c(1, 4, 5))
layout(m)

par(mar = c(3, 4, 1, 1), oma = c(0, 0, 1, 0), mgp=c(3,1,0))

highlight <- subset(mydata, pH_selection_pval < cut_off)
manhattan(mydata, chr="CHR", bp="POS", p="pH_selection_pval",
    genomewideline = FALSE,
    suggestiveline = FALSE, 
    highlight= highlight,
    ylim=c(0, (-log10(min(mydata$pH_selection_pval, na.rm=TRUE)))*1.1),
    xlab="", cex.axis=0.5, 
    ylab="", xaxt="n", yaxt="n")

#axis(1, mgp=c(1.8, .2, 0), cex.axis=0.6) # second is tick mark labels
title(ylab=expression(-log[10](italic(p))), line=2, cex.lab=0.8)
title(xlab="Scaffold", line=0.5, cex.lab=0.8)
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7)

box()

abline(h=-log10(cut_off), lty=2, col="firebrick3", lwd=1.5)

# add label
mtext(text="A",
        side=3, line=0,
             cex=1.5, 
            at=par("usr")[1]-0.06*diff(par("usr")[1:2]), outer=FALSE)

## plot 2 p value
#d <- d2
d <- subset(d, pos > 34300 & pos < 34600)

#par(mfrow = c(2, 2))
par(mar=c(3, 4, 1, 1), mgp=c(3, 1, 0), las=0)

#pval plot:
plot(d$pos, -log10(d$logp), xlim=c(34350,34510), pch=19, col="black",
    xlab="",
    ylab="",xaxt="n",yaxt="n", cex=0.8)
#points(x=highlight$pos, y=-log10(highlight$logp), pch=19, cex=1.1, col="green")
axis(1, mgp=c(1.8, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7)
abline(h=-log10(cut_off), lty=2, col="firebrick3", lwd=1.5)
title(ylab=expression(-log[10](italic(p))), line=2, cex.lab=0.8)

# add label
mtext(text="B",
        side=3, line=0,
             cex=1.5, 
            at=par("usr")[1]-0.13*diff(par("usr")[1:2]))


### af change

par(mar=c(4, 4, 0, 1), mgp=c(3, 1, 0), las=0)

# af change plot:
plot(d$pos, d$D7_7_delta, xlim=c(34350,34510), ylim=c(0,0.22), pch=19, col="firebrick3",
    xlab="",
    ylab="",xaxt="n",yaxt="n")
axis(1, mgp=c(1.8, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7)
#title(ylab=expression(atop(paste(Delta, " allele frequency"), "vs. D1 control")), line=1.9, cex.lab=0.8)
mtext(side=2,expression(paste(Delta, " allele frequency")), line=1.9, cex=0.5)
mtext(side=2, "vs. D1 control", line=1.2, cex=0.5)

title(xlab="Scaffold542: position in BP", line=1.9, cex.lab=0.8)

lines((d$pos), d$D7_7_delta, col="firebrick3")
points((d$pos), d$D1_7_delta, col="springgreen4", pch=17)
lines((d$pos), d$D1_7_delta, col="springgreen4")
points((d$pos), d$D7_8_delta, col="royalblue3", pch=15)
lines((d$pos), d$D7_8_delta, col="royalblue3")

mtext(text="C",
        side=3, line=0,
             cex=1.5, 
            at=par("usr")[1]-0.13*diff(par("usr")[1:2]))

######## plot 3

# currently, this is 
new.dat <- subset(mydata, CHR > 0 & CHR < 2)

d=data.frame(CHR=new.dat$CHR, BP=new.dat$POS, P=new.dat$pH_selection_pval, SNP=new.dat$SNP,
            D7_7_delta = new.dat$D7_7_delta,D1_7_delta = new.dat$D1_7_delta, D7_8_delta=new.dat$D7_8_delta) 
    
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

d <- subset(d, pos > 665365 & pos < 665510)


### second set of plots
#par(mfrow = c(2, 2))
par(mar=c(3, 4, 1, 1), mgp=c(3, 1, 0), las=0)

#pval plot:
plot(d$pos, -log10(d$logp), xlim=c(665365,665510), pch=19, col="black",
    xlab="",
    ylab="",xaxt="n",yaxt="n", cex=0.8)
#points(x=highlight$pos, y=-log10(highlight$logp), pch=19, cex=1.1, col="green")
axis(1, mgp=c(1.8, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7)
abline(h=-log10(cut_off), lty=2, col="firebrick3", lwd=1.5)
title(ylab=expression(-log[10](italic(p))), line=2, cex.lab=0.8)


mtext(text="D",
        side=3, line=0,
             cex=1.5, 
            at=par("usr")[1]-0.13*diff(par("usr")[1:2]))

# af change plot:
par(mar=c(4, 4, 0, 1), mgp=c(3, 1, 0), las=0)

plot(d$pos, d$D7_7_delta, xlim=c(665365,665510), ylim=c(0,0.22), pch=19, col="firebrick3",
    xlab="",
    ylab="",xaxt="n",yaxt="n")
axis(1, mgp=c(1.8, .5, 0), cex.axis=0.7) # second is tick mark labels
axis(2, mgp=c(1.8, .4, 0), cex.axis=0.7)
mtext(side=2,expression(paste(Delta, " allele frequency")), line=1.9, cex=0.5)
mtext(side=2, "vs. D1 control", line=1.2, cex=0.5)
title(xlab="Scaffold1: position in BP", line=1.9, cex.lab=0.8)

lines((d$pos), d$D7_7_delta, col="firebrick3")
points((d$pos), d$D1_7_delta, col="springgreen4", pch=17)
lines((d$pos), d$D1_7_delta, col="springgreen4")
points((d$pos), d$D7_8_delta, col="royalblue3", pch=15)
lines((d$pos), d$D7_8_delta, col="royalblue3")

mtext(text="E",
        side=3, line=0,
             cex=1.5, 
            at=par("usr")[1]-0.13*diff(par("usr")[1:2]))

# add legend. basically overlaying new empty plot
par(fig = c(0.05, 1, 0, 0.11), oma = c(0, 0, .15, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0,0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("center", legend= c("D1 pH 7.5", "D7 pH 8.0", "D7 pH 7.5"), xpd = TRUE, 
    horiz = FALSE, inset = c(0, 0), 
    bty="n",pch = c(15,17,19), col = c("springgreen4", "royalblue3", "firebrick3"), lwd=1, cex=0.8, pt.cex=1.2)

dev.off()

