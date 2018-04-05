#library(qqman)

mydata <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
mydata$CHR <- as.numeric(gsub("Scaffold", "", mydata$CHROM))
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

cut_off <- 0.001

#allele frequencies to plot panel
#calc stderror of mean af estimates

gps <- c("D1_7", "D7_7", "D7_8")

mydata$D7_7_delta <- abs(mydata$D1_8_mean-mydata$D7_7_mean)
mydata$D7_8_delta <- abs(mydata$D1_8_mean-mydata$D7_8_mean)

new.dat <- mydata

d=data.frame(CHR=new.dat$CHR, BP=new.dat$POS, P=-log10(new.dat$pH7_selection_qval), 
            SNP=new.dat$SNP)

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
    ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
    labs <- unique(d$CHR)
    }

    # Initialize plot
    xmax = ceiling(max(d$pos) * 1.005)
    xmin = floor(max(d$pos) * -0.005)

# define plot parameters:
    def_args <- list(xaxt='n', yaxt='n', xaxs='i', yaxs='i', las=1, pch=20,
                     xlim=c(xmin,xmax), ylim=c(0,24.42359),
                     xlab="", ylab="", cex=0.9, bty="n")



png("~/urchin_af/figures/Fig_S2_manhattan.png",width=6.65, height= 4, units="in", res=300)

m <- cbind(c(1,2), c(1, 2))
layout(m)

par(mar = c(3, 4, 1, 1), oma = c(0, 0, 1, 0), mgp=c(3,1,0))

selected_7 <- mydata[which(mydata$pH7_selection_qval < cut_off & mydata$pH8_selection_qval >= cut_off),]
selected_8 <- mydata[which(mydata$pH8_selection_qval < cut_off & mydata$pH7_selection_qval >= cut_off),]
selected_both <- mydata[which(mydata$pH8_selection_qval < cut_off & mydata$pH7_selection_qval < cut_off),]

    do.call("plot", c(NA, def_args))    
    # Add an axis. 
    labs[1:length(labs)] <- ""
    #axis(1, at=ticks, labels=rep("", length(ticks)), cex.axis=0.8, mgp=c(0,0.4,0.1))
    #axis(1, at=ticks[(seq(2,24,2))], labels=labs[(seq(2,24,2))], cex.axis=0.8, mgp=c(0,0.4,0.1))
    my_labs <- seq(0,24.42359 , by = 4)
    axis(side = 2, at = my_labs, labels = my_labs, cex.axis=0.7,mgp=c(0,0.4,0.1))


    # Create a vector of alternatiting colors    
    col=c("gray10", "gray60")
    col = rep_len(col, max(d$index)) 

    # Add points to the plot
    if (nchr==1) {
        with(d, points(pos, logp, pch=20, col=col[1]))
    } else {
        icol=1
        for (i in unique(d$index)) {
        points(d[d$index==i,"pos"], d[d$index==i,"logp"], col=col[icol], pch=20)
            icol=icol+1
        }
    }

    # Highlight snps from a character vector
    d.highlight <- d[which(d$SNP %in% selected_7$SNP),]
    with(d.highlight, points(pos, logp, col="black", bg="firebrick1", pch=21)) 

    d.highlight <- d[which(d$SNP %in% selected_both$SNP),]
    with(d.highlight, points(pos, logp, col="black", bg="darkorchid1", pch=24, cex=1)) 

title(ylab=expression(-log[10](italic(p))), line=1.5, cex.lab=0.8)
title(xlab="Scaffold", line=0.5, cex.lab=0.8)

# add label
mtext(text="A",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.06*diff(par("usr")[1:2]), outer=FALSE)


########
## pH 8.0 plot
########

d=data.frame(CHR=new.dat$CHR, BP=new.dat$POS, P=-log10(new.dat$pH8_selection_qval), 
            SNP=new.dat$SNP)

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
    xmax = ceiling(max(d$pos) * 1.005)
    xmin = floor(max(d$pos) * -0.005)

# define plot parameters:
    def_args <- list(xaxt='n', yaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                     xlim=c(xmin,xmax), ylim=c(0,24.42359),
                     xlab="", ylab="")

#par(mar = c(3, 4, 1, 1), oma = c(0, 0, 1, 0), mgp=c(3,1,0))

    do.call("plot", c(NA, def_args))    
    # Add an axis. 
    labs[1:length(labs)] <- ""
   # axis(1, at=ticks, labels=rep("", length(ticks)), cex.axis=0.8, mgp=c(0,0.4,0.1))
    #axis(1, at=ticks[(seq(2,24,2))], labels=labs[(seq(2,24,2))], cex.axis=0.8, mgp=c(0,0.4,0.1))
    my_labs <- seq(0,24.42359, by = 4)
    axis(side = 2, at = my_labs, labels = my_labs, cex.axis=0.7,mgp=c(0,0.4,0.1))


    # Create a vector of alternatiting colors    
    col=c("gray10", "gray60")
    col = rep_len(col, max(d$index)) 

    # Add points to the plot
    if (nchr==1) {
        with(d, points(pos, logp, pch=20, col=col[1]))
    } else {
        icol=1
        for (i in unique(d$index)) {
        points(d[d$index==i,"pos"], d[d$index==i,"logp"], col=col[icol], pch=20)
            icol=icol+1
        }
    }

    # Highlight snps from a character vector
    d.highlight <- d[which(d$SNP %in% selected_8$SNP),]
    with(d.highlight, points(pos, logp, col="black", bg="royalblue1", pch=22, cex=1)) 

    d.highlight <- d[which(d$SNP %in% selected_both$SNP),]
    with(d.highlight, points(pos, logp, col="black", bg="darkorchid1", pch=24, cex=1)) 

title(ylab=expression(-log[10](italic(p))), line=1.5, cex.lab=0.8)
title(xlab="Scaffold", line=0.5, cex.lab=0.8)

# add label
mtext(text="B",
        side=3, line=0,
             cex=1.5,
            at=par("usr")[1]-0.06*diff(par("usr")[1:2]), outer=FALSE)

par(fig = c(0.00, 0.97, 0, 0.55), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0,0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("topright", legend= c("pH 7.5", "pH 8.0", "overlapping variants"),
    horiz = FALSE,
    pch = c(21,22,24), col=c("black", "black", "black"), pt.bg = c( "firebrick1", "royalblue1", "darkorchid1"), 
    cex=0.9, pt.cex=1.3)



dev.off()
