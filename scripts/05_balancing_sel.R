brary(stringr)
library(ggplot2)
library(gridExtra)
library(MASS)

dat <- read.table("~/urchin_af/variants/urchin_ann.vcf", stringsAsFactors=FALSE)

# loop over entire data set. count up annotattions. Make sure I'm pulling out what I think I am
eff <- strsplit(as.character(dat$V8), split=",", fixed=TRUE)

one=0
two=0
three=0
four=0
five=0
six=0
seven=0
eight=0

for (i in 1:length(eff)){
    if(length(eff[[i]]) == 1){
        one=one+1
    }
    if(length(eff[[i]]) == 2){
        two=two+1
    }
    if(length(eff[[i]]) == 3){
        three=three+1
    }
    if(length(eff[[i]]) == 4){
        four=four+1
    }
    if(length(eff[[i]]) == 5){
        five=five+1
    }
    if(length(eff[[i]]) == 6){
        six=six+1
    }
    if(length(eff[[i]]) == 7){
        seven=seven+1
    }
    if(length(eff[[i]]) == 8){
        eight=eight+1
    }
}

#sanity check
(one+two+three+four+five+six+seven+eight) == nrow(dat)


# it is possible that snps receive multiple annotations
# I think in this case, it is best to take the most "severe"
# order is: HIGH MODERATE LOW MODIFIER

## from snpeff manual:
# Effect sort order. When multiple effects are reported, SnpEff sorts the effects the following way:
    # Putative impact: Effects having higher putative impact are first.
    # Effect type: Effects assumed to be more deleterious effects first.
    # Canonical trancript before non-canonical.
    # Marker genomic coordinates (e.g. genes starting before first).

# bc of this, should be safe to simply take the first annotation,
    # ie, the one with the highest potential impact

# make dataframe
out <- dat

for (i in 1:length(eff)){
    if(length(eff[[i]]) == 1){
        out$V8[i] <- eff[[i]]
    }
    else{
        # pull out the effect of each variant
        # they're ordered by severity
        out$V8[i] <- eff[[i]][1]
        }
    }

# add snp id
out$SNP <- paste(out$V1, out$V2, sep=":")

out_ann <- data.frame(SNP=out$SNP, ANN=out$V8)
# split annotation
ann_split <- (str_split_fixed(out$V8, "\\|", n=16))
out_ann <-(cbind(out$SNP, ann_split))

colnames(out_ann) <- c("SNP","Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID", "Feature_Type","Feature_ID", "Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos",  "cDNA.length","CDS.pos","Distance","WARNINGS")

########################################
# assign functional categories
########################################

#intergenic: downstream_gene_variant, intergenic_region, upstream_gene_variant
#intron: intron_variant, splice_region_variant&intron_variant,
    # splice_acceptor_variant&intron_variant,
    # splice_donor_variant&splice_region_variant&intron_variant,
    # splice_region_variant, splice_region_variant&intron_variant,
    # splice_region_variant&synonymous_variant
#synonymous: synonymous_variant, stop_retained_variant
#non-synonymous: missense_variant, stop_gained, stop_lost,
    # missense_variant&splice_region_variant, start_lost, stop_gained&splice_region_variant,
    # stop_lost&splice_region_variant, initiator_codon_variant
#3' UTR: 3_prime_UTR_variant
# unclear: non_coding_transcript_exon_variant,
    # splice_region_variant&non_coding_transcript_exon_variant
    # splice_region_variant&stop_retained_variant

# looking at column Feature_Type

# read in cmh results so can identify sig changes:

cmh <- read.table("~/urchin_af/analysis/cmh.out.txt", header=TRUE)
cmh$SNP <- paste(cmh$CHROM, cmh$POS, sep=":")
new <- as.data.frame(out_ann)

# assign significance to new
new$sig_pH7 <- cmh$pH7_sig[match(new$SNP, cmh$SNP)]
new$sig_pH8 <- cmh$pH8_sig[match(new$SNP, cmh$SNP)]

# assign class  ie, intergenic, etc
new$class <- c(NA)

new$class[which(new$Annotation == "downstream_gene_variant" |
    new$Annotation == "intergenic_region" |
    new$Annotation == "upstream_gene_variant")] <- c("intergenic")
new$class[which(new$Annotation == "intron_variant" |
    new$Annotation == "splice_region_variant&intron_variant"|
    new$Annotation == "splice_acceptor_variant&intron_variant"|
    new$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"|
    new$Annotation == "splice_region_variant"|
    new$Annotation == "splice_region_variant&intron_variant" |
    new$Annotation == "splice_region_variant&synonymous_variant"|
    new$Annotation == "splice_region_variant&non_coding_transcript_exon_variant"|
    new$Annotation == "non_coding_transcript_exon_variant"|
    new$Annotation == "splice_donor_variant&intron_variant")] <- c("intron")
new$class[which(new$Annotation == "synonymous_variant" |
    new$Annotation == "stop_retained_variant"|
    new$Annotation == "splice_region_variant&stop_retained_variant")] <- c("synonymous")
new$class[which(new$Annotation == "missense_variant" |
    new$Annotation == "stop_gained" |
    new$Annotation == "stop_lost"|
    new$Annotation == "missense_variant&splice_region_variant"|
    new$Annotation == "start_lost"|
    new$Annotation == "stop_gained&splice_region_variant"|
    new$Annotation == "stop_lost&splice_region_variant"|
    new$Annotation == "initiator_codon_variant"|
    new$Annotation == "3_prime_UTR_variant")] <- c("non-synonymous")

ens <- as.data.frame(matrix(ncol=3, nrow=4))
colnames(ens) <- c("gff", "SNP_class", "count")
ens$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous")
ens$gff <- "ensembl"
ens$count[1] <- length(which(new$class == "intergenic"))
ens$count[2] <- length(which(new$class == "intron"))
ens$count[3] <- length(which(new$class == "synonymous"))
ens$count[4] <- length(which(new$class == "non-synonymous"))

ens$proportion <- ens$count/sum(ens$count)

dat_count <- ens

pdf("~/urchin_af/figures/snp_class.pdf", height=7, width=7)

a <- ggplot(data=dat_count, aes(x=SNP_class, y=proportion, group="SNP_class")) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('gray50'))
a

dev.off()

#############
# split into selected and non
#############

###############
# pH 7.5 selected loci
ens.sel_7 <- as.data.frame(matrix(ncol=3, nrow=4))
colnames(ens.sel_7) <- c("gff", "SNP_class", "count")
ens.sel_7$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous")
ens.sel_7$gff <- "selected_pH7"
new.sel_7 <- new[which(new$sig_pH7 ==TRUE),]

ens.sel_7$count[1] <- length(which(new.sel_7$class == "intergenic"))
ens.sel_7$count[2] <- length(which(new.sel_7$class == "intron"))
ens.sel_7$count[3] <- length(which(new.sel_7$class == "synonymous"))
ens.sel_7$count[4] <- length(which(new.sel_7$class == "non-synonymous"))

ens.sel_7$proportion <- ens.sel_7$count/sum(ens.sel_7$count)

###############
# pH 8.0 selected loci
ens.sel_8 <- as.data.frame(matrix(ncol=3, nrow=4))
colnames(ens.sel_8) <- c("gff", "SNP_class", "count")
ens.sel_8$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous")
ens.sel_8$gff <- "selected_pH8"
new.sel_8 <- new[which(new$sig_pH8 ==TRUE),]

ens.sel_8$count[1] <- length(which(new.sel_8$class == "intergenic"))
ens.sel_8$count[2] <- length(which(new.sel_8$class == "intron"))
ens.sel_8$count[3] <- length(which(new.sel_8$class == "synonymous"))
ens.sel_8$count[4] <- length(which(new.sel_8$class == "non-synonymous"))

ens.sel_8$proportion <- ens.sel_8$count/sum(ens.sel_8$count)

###############
# neutral loci
ens.neut <- as.data.frame(matrix(ncol=3, nrow=4))
colnames(ens.neut) <- c("gff", "SNP_class", "count")
ens.neut$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous")
ens.neut$gff <- "neutral"
new.neut <- new[which(new$sig_pH8 ==FALSE & new$sig_pH7 ==FALSE),]

ens.neut$count[1] <- length(which(new.neut$class == "intergenic"))
ens.neut$count[2] <- length(which(new.neut$class == "intron"))
ens.neut$count[3] <- length(which(new.neut$class == "synonymous"))
ens.neut$count[4] <- length(which(new.neut$class == "non-synonymous"))

ens.neut$proportion <- ens.neut$count/sum(ens.neut$count)

ens_all <- rbind(ens.neut, ens.sel_8, ens.sel_7)

a <- ggplot(data=ens_all, aes(x=SNP_class, y=proportion, fill=gff)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('azure4','tomato3','royalblue4')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.title=element_blank())
a


###################################################
###################################################
##################################################
# write output
###################################################
###################################################
###################################################

write.table(file="~/urchin_af/analysis/cmh.annotations.out", new,
    col.name=TRUE, quote=FALSE, row.name=FALSE, sep= "\t" )

###################################################
###################################################
###################################################
### I would be super cool to see (with either/both) the
# distribution of NS across MAF bins for selected and neutral loci.
###################################################
###################################################
###################################################

# need to:
    # identify the Non-syn variants
    # pull out selected vs neutral
    # bin into maf
        # both folded and not
    # repeat for transcriptome and ensembl


# read in allele frequency data
mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE, stringsAsFactors=FALSE)
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

new.sel_7 <- new[which(new$sig_pH7 ==TRUE),]
new.sel_8 <- new[which(new$sig_pH8 ==TRUE),]
new.neut <- new[which(new$sig_pH8 ==FALSE & new$sig_pH7 ==FALSE),]


# fold af
mydata$folded_af = unlist(lapply(mydata$D1_8_mean,function(x)
          ifelse(x > 0.5, (1-x), x)), use.names=FALSE)


# ns variants
ns.sel_7 <- new.sel_7[which(new.sel_7$class == "non-synonymous"),]
ns.sel_8 <- new.sel_8[which(new.sel_8$class == "non-synonymous"),]
ns.neut <- new.neut[which(new.neut$class == "non-synonymous"),]

# syn variants
syn.sel_7 <- new.sel_7[which(new.sel_7$class == "synonymous"),]
syn.sel_8 <- new.sel_8[which(new.sel_8$class == "synonymous"),]
syn.neut <- new.neut[which(new.neut$class == "synonymous"),]

# introns
intron.sel_7 <- new.sel_7[which(new.sel_7$class == "intron"),]
intron.sel_8 <- new.sel_8[which(new.sel_8$class == "intron"),]
intron.neut <- new.neut[which(new.neut$class == "intron"),]

#intergenic
intergen.sel_7 <- new.sel_7[which(new.sel_7$class == "intergenic"),]
intergen.sel_8 <- new.sel_8[which(new.sel_8$class == "intergenic"),]
intergen.neut <- new.neut[which(new.neut$class == "intergenic"),]


# unfolded allele freqs are col: af_out
# merge datasets:
ns.sel_7.merge <- merge(mydata,ns.sel_7, by="SNP" )
ns.sel_8.merge <- merge(mydata,ns.sel_8, by="SNP" )
ns.neut.merge <- merge(mydata,ns.neut, by="SNP" )
ns.all <- rbind(ns.neut.merge,ns.sel_7.merge,ns.sel_8.merge)

syn.sel_7.merge <- merge(mydata,syn.sel_7, by="SNP" )
syn.sel_8.merge <- merge(mydata,syn.sel_8, by="SNP" )
syn.neut.merge <- merge(mydata,syn.neut, by="SNP" )
syn.all <- rbind(syn.neut.merge,syn.sel_7.merge,syn.sel_8.merge)

intron.sel_7.merge <- merge(mydata,intron.sel_7, by="SNP" )
intron.sel_8.merge <- merge(mydata,intron.sel_8, by="SNP" )
intron.neut.merge <- merge(mydata,intron.neut, by="SNP" )
intron.all <- rbind(intron.neut.merge,intron.sel_7.merge,intron.sel_8.merge)

intergen.sel_7.merge <- merge(mydata,intergen.sel_7, by="SNP" )
intergen.sel_8.merge <- merge(mydata,intergen.sel_8, by="SNP" )
intergen.neut.merge <- merge(mydata,intergen.neut, by="SNP" )
intergen.all <- rbind(intergen.neut.merge, intergen.sel_7.merge,intergen.sel_8.merge)


# unfolded
ns.neut.merge$unfold_bin <- cut(ns.neut.merge$af_out, breaks=seq(0,1, 0.05))
ns.sel_7.merge$unfold_bin <- cut(ns.sel_7.merge$af_out, breaks=seq(0,1, 0.05))
ns.sel_8.merge$unfold_bin <- cut(ns.sel_8.merge$af_out, breaks=seq(0,1, 0.05))
ns.all$unfold_bin <- cut(ns.all$af_out, breaks=seq(0,1, 0.05))

syn.neut.merge$unfold_bin <- cut(syn.neut.merge$af_out, breaks=seq(0,1, 0.05))
syn.sel_7.merge$unfold_bin <- cut(syn.sel_7.merge$af_out, breaks=seq(0,1, 0.05))
syn.sel_8.merge$unfold_bin <- cut(syn.sel_8.merge$af_out, breaks=seq(0,1, 0.05))
syn.all$unfold_bin <- cut(syn.all$af_out, breaks=seq(0,1, 0.05))

intron.neut.merge$unfold_bin <- cut(intron.neut.merge$af_out, breaks=seq(0,1, 0.05))
intron.sel_7.merge$unfold_bin <- cut(intron.sel_7.merge$af_out, breaks=seq(0,1, 0.05))
intron.sel_8.merge$unfold_bin <- cut(intron.sel_8.merge$af_out, breaks=seq(0,1, 0.05))
intron.all$unfold_bin <- cut(intron.all$af_out, breaks=seq(0,1, 0.05))

intergen.neut.merge$unfold_bin <- cut(intergen.neut.merge$af_out, breaks=seq(0,1, 0.05))
intergen.sel_7.merge$unfold_bin <- cut(intergen.sel_7.merge$af_out, breaks=seq(0,1, 0.05))
intergen.sel_8.merge$unfold_bin <- cut(intergen.sel_8.merge$af_out, breaks=seq(0,1, 0.05))
intergen.all$unfold_bin <- cut(intergen.all$af_out, breaks=seq(0,1, 0.05))

####### Unfolded #########

bin <- c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25",
    "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5",
    "0.5-0.55", "0.55-0.6", "0.6-0.65", "0.65-0.7", "0.7-0.75",
    "0.75-0.8", "0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-1")

# make barplot of neutral vs selected
ns_neut_unfoldfreq <- table(ns.neut.merge$unfold_bin)/sum(table(ns.neut.merge$unfold_bin))
ns_sel_7_unfoldfreq <- table(ns.sel_7.merge$unfold_bin)/sum(table(ns.sel_7.merge$unfold_bin))
ns_sel_8_unfoldfreq <- table(ns.sel_8.merge$unfold_bin)/sum(table(ns.sel_8.merge$unfold_bin))
ns_all_unfoldfreq <- table(ns.all$unfold_bin)/sum(table(ns.all$unfold_bin))

syn_neut_unfoldfreq <- table(syn.neut.merge$unfold_bin)/sum(table(syn.neut.merge$unfold_bin))
syn_sel_7_unfoldfreq <- table(syn.sel_7.merge$unfold_bin)/sum(table(syn.sel_7.merge$unfold_bin))
syn_sel_8_unfoldfreq <- table(syn.sel_8.merge$unfold_bin)/sum(table(syn.sel_8.merge$unfold_bin))
syn_all_unfoldfreq <- table(syn.all$unfold_bin)/sum(table(syn.all$unfold_bin))

intron_neut_unfoldfreq <- table(intron.neut.merge$unfold_bin)/sum(table(intron.neut.merge$unfold_bin))
intron_sel_7_unfoldfreq <- table(intron.sel_7.merge$unfold_bin)/sum(table(intron.sel_7.merge$unfold_bin))
intron_sel_8_unfoldfreq <- table(intron.sel_8.merge$unfold_bin)/sum(table(intron.sel_8.merge$unfold_bin))
intron_all_unfoldfreq <- table(intron.all$unfold_bin)/sum(table(intron.all$unfold_bin))

intergen_neut_unfoldfreq <- table(intergen.neut.merge$unfold_bin)/sum(table(intergen.neut.merge$unfold_bin))
intergen_sel_7_unfoldfreq <- table(intergen.sel_7.merge$unfold_bin)/sum(table(intergen.sel_7.merge$unfold_bin))
intergen_sel_8_unfoldfreq <- table(intergen.sel_8.merge$unfold_bin)/sum(table(intergen.sel_8.merge$unfold_bin))
intergen_all_unfoldfreq <- table(intergen.all$unfold_bin)/sum(table(intergen.all$unfold_bin))

# merge into tables for plotting
ns.sel_7.dat <- data.frame(class=rep("pH 7.5 selected", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_sel_7_unfoldfreq))
ns.sel_8.dat <- data.frame(class=rep("pH 8.0 selected", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_sel_8_unfoldfreq))
ns.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_neut_unfoldfreq))
ns.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_all_unfoldfreq))
ns.fq.dat <- rbind(ns.sel_7.dat,ns.sel_8.dat, ns.neut.dat, ns.all.dat)

syn.sel_7.dat <- data.frame(class=rep("pH 7.5 selected", length(bin)), bin=bin,
            allele_frequency=as.vector(syn_sel_7_unfoldfreq))
syn.sel_8.dat <- data.frame(class=rep("pH 8.0 selected", length(bin)), bin=bin,
            allele_frequency=as.vector(syn_sel_8_unfoldfreq))
syn.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(syn_neut_unfoldfreq))
syn.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(syn_all_unfoldfreq))
syn.fq.dat <- rbind(syn.sel_7.dat,syn.sel_8.dat, syn.neut.dat, syn.all.dat)

intron.sel_7.dat <- data.frame(class=rep("pH 7.5 selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intron_sel_7_unfoldfreq))
intron.sel_8.dat <- data.frame(class=rep("pH 8.0 selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intron_sel_8_unfoldfreq))
intron.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(intron_neut_unfoldfreq))
intron.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(intron_all_unfoldfreq))
intron.fq.dat <- rbind(intron.sel_7.dat,intron.sel_8.dat, intron.neut.dat, intron.all.dat)

intergen.sel_7.dat <- data.frame(class=rep("pH 7.5 selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intergen_sel_7_unfoldfreq))
intergen.sel_8.dat <- data.frame(class=rep("pH 8.0 selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intergen_sel_8_unfoldfreq))
intergen.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(intergen_neut_unfoldfreq))
intergen.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(intergen_all_unfoldfreq))
intergen.fq.dat <- rbind(intergen.sel_7.dat,intergen.sel_8.dat, intergen.neut.dat, intergen.all.dat)



a <- ggplot(data=ns.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'darkolivegreen3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 0.45) +
  guides(size = FALSE) +
  ggtitle("non-synonymous")+
  ylab("relative frequency") +
  theme(legend.title=element_blank())

b <- ggplot(data=syn.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'darkolivegreen3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 0.45) +
  guides(size = FALSE) +
  ggtitle("synonymous")+
  ylab("relative frequency") +
  theme(legend.title=element_blank())

c <- ggplot(data=intron.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'darkolivegreen3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 0.45) +
  guides(size = FALSE) +
  ggtitle("intron")+
  ylab("relative frequency") +
  theme(legend.title=element_blank())

d <- ggplot(data=intergen.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'darkolivegreen3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("intergenic")+
  ylim(0, 0.45) +
  guides(size = FALSE) +
  ylab("relative frequency") +
  theme(legend.title=element_blank())

library(ggpubr)

png("~/urchin_af/figures/AF_categories_unfolded_bin.png", height=7, res=300, units="in", width=10)
ggarrange(a, b, c, d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()


#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(a)

p3 <- grid.arrange(arrangeGrob(a + theme(legend.position="none"),
                         b + theme(legend.position="none"),
                         c + theme(legend.position="none"),
                         d + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=2,heights=c(10, 1))

png("~/urchin_af/figures/AF_categories_unfolded_bin.png", height=7, res=300, units="in", width=10)
p3
dev.off()

### Some stats

ks.test(ns.sel_7.merge$af_out, ns.sel_8.merge$af_out)
ks.test(syn.sel_7.merge$af_out, syn.sel_8.merge$af_out)
ks.test(intron.sel_7.merge$af_out, intron.sel_8.merge$af_out)
ks.test(intergen.sel_7.merge$af_out, intergen.sel_8.merge$af_out)

