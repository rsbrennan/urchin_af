d1_8 <- read.table("~/urchin_af/analysis/D1_8.pop.windowed.pi", header=TRUE)
d1_7 <- read.table("~/urchin_af/analysis/D1_7.pop.windowed.pi", header=TRUE)
d7_7 <- read.table("~/urchin_af/analysis/D7_7.pop.windowed.pi", header=TRUE)
d7_8 <- read.table("~/urchin_af/analysis/D7_8.pop.windowed.pi", header=TRUE)

d7_7.sel <- read.table("~/urchin_af/analysis/D7_7.pop.selected.windowed.pi", header=TRUE)
d1_8.sel <- read.table("~/urchin_af/analysis/D1_8.pop.selected.windowed.pi", header=TRUE)
d7_8.sel <- read.table("~/urchin_af/analysis/D7_8.pop.selected.windowed.pi", header=TRUE)

d1_8$samp <- c("d1_8")
d1_7$samp <- c("d1_7")
d7_7$samp <- c("d7_7")
d7_8$samp <- c("d7_8")
d1_8.sel$samp <- c("d1_8.sel")
d7_7.sel$samp <- c("d7_7.sel")

dat <- rbind(d1_7, d1_8, d7_7, d7_8)

test <- d1_8[sample(1:nrow(d1_8),nrow(d1_8.sel)),]
dat <- rbind(d1_8.sel, test)
boxplot(dat$PI ~ dat$samp)

hist(d1_8.sel$PI, breaks=100, col="grey", freq=FALSE)
hist(d7_7.sel$PI, breaks=100, col="red", freq=FALSE, add=T)

mean(d1_8$PI)
mean(d1_7$PI)
mean(d7_7$PI)
mean(d7_8$PI)

# look at windows where n > 1
mean(d1_8$PI[which(d1_8$N_VARIANTS >1)])
mean(d7_7$PI[which(d7_7$N_VARIANTS >1)])
mean(d1_7$PI[which(d1_7$N_VARIANTS >1)])
mean(d7_8$PI[which(d7_8$N_VARIANTS >1)])

ks.test(x=d7_8$PI[which(d7_8$N_VARIANTS >1)], y=d7_7$PI[which(d7_7$N_VARIANTS >1)])
ks.test(x=d1_8$PI[which(d1_8$N_VARIANTS >1)], y=d1_7$PI[which(d1_7$N_VARIANTS >1)])
ks.test(x=d1_8$PI, y=d7_7$PI)
ks.test(x=d1_8$PI, y=d7_8$PI)


dat.new <- dat[which(dat$N_VARIANTS > 1),]
boxplot(dat.new$PI ~ dat.new$samp)

ks.test(x=d7_8$PI, y=d7_7$PI)
ks.test(x=d1_8$PI, y=d1_7$PI)
ks.test(x=d1_8$PI, y=d7_7$PI)
ks.test(x=d1_8$PI, y=d7_8$PI)
