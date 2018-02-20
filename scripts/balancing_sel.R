# snp assignments

library(stringr)
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
        #change <- lapply(strsplit(eff[[i]], split="|", fixed=TRUE), '[[', 2)
        #change.eff <- lapply(strsplit(eff[[i]], split="|", fixed=TRUE), '[[', 3)
        #check if they all match. if so, just take the first one
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

new <- as.data.frame(out_ann)
new$sig <- cmh$pH_sig
new$class <- c(NA) # also want to add in "class" ie, intergenic, etc

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
# selected loci
ens.sel <- as.data.frame(matrix(ncol=3, nrow=4))
colnames(ens.sel) <- c("gff", "SNP_class", "count")
ens.sel$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous")
ens.sel$gff <- "selected"
new.sel <- new[which(new$sig ==TRUE),]

ens.sel$count[1] <- length(which(new.sel$class == "intergenic"))
ens.sel$count[2] <- length(which(new.sel$class == "intron"))
ens.sel$count[3] <- length(which(new.sel$class == "synonymous"))
ens.sel$count[4] <- length(which(new.sel$class == "non-synonymous"))

ens.sel$proportion <- ens.sel$count/sum(ens.sel$count)

###############
# neutral loci
ens.neut <- as.data.frame(matrix(ncol=3, nrow=4))
colnames(ens.neut) <- c("gff", "SNP_class", "count")
ens.neut$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous")
ens.neut$gff <- "neutral"
new.neut <- new[which(new$sig ==FALSE),]

ens.neut$count[1] <- length(which(new.neut$class == "intergenic"))
ens.neut$count[2] <- length(which(new.neut$class == "intron"))
ens.neut$count[3] <- length(which(new.neut$class == "synonymous"))
ens.neut$count[4] <- length(which(new.neut$class == "non-synonymous"))

ens.neut$proportion <- ens.neut$count/sum(ens.neut$count)

ens_all <- rbind(ens.neut, ens.sel)

a <- ggplot(data=ens_all, aes(x=SNP_class, y=proportion, fill=gff)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.title=element_blank())
a

# make contingency table

r1 <- c(4375, 40709-4375)
r2 <- c(98, 785-98)

din <- matrix(c(r1, r2),
                nrow=2,
                byrow=TRUE)
rownames(din) <- c("neutral", "selected")
colnames(din) <- c("non-synonymous", "all")

chisq.test(din,correct=FALSE)



####### 
# write output

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
mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")
new.sel <- new[which(new$sig ==TRUE),]
new.neut <- new[which(new$sig ==FALSE),]

# fold af
mydata$folded_af = unlist(lapply(mydata$D1_8_mean,function(x)
          ifelse(x > 0.5, (1-x), x)), use.names=FALSE)


# ns variants
ns.sel <- new.sel[which(new.sel$class == "non-synonymous",]
ns.neut <- new.neut[which(new.neut$class == "non-synonymous",]

# syn variants
syn.sel <- new.sel[new.sel$class == "synonymous",]
syn.neut <- new.neut[new.neut$class == "synonymous",]

# introns
intron.sel <- new.sel[new.sel$class == "intron",]
intron.neut <- new.neut[new.neut$class == "intron",]

#intergenic
intergen.sel <- new.sel[new.sel$class == "intergenic",]
intergen.neut <- new.neut[new.neut$class == "intergenic",]


# unfolded allele freqs are col: af_out
# merge datasets:
ns.sel.merge <- merge(mydata,ns.sel, by="SNP" )
ns.neut.merge <- merge(mydata,ns.neut, by="SNP" )
ns.all <- rbind(ns.neut.merge,ns.sel.merge)

syn.sel.merge <- merge(mydata,syn.sel, by="SNP" )
syn.neut.merge <- merge(mydata,syn.neut, by="SNP" )
syn.all <- rbind(syn.neut.merge,syn.sel.merge)

intron.sel.merge <- merge(mydata,intron.sel, by="SNP" )
intron.neut.merge <- merge(mydata,intron.neut, by="SNP" )
intron.all <- rbind(intron.neut.merge,intron.sel.merge)

intergen.sel.merge <- merge(mydata,intergen.sel, by="SNP" )
intergen.neut.merge <- merge(mydata,intergen.neut, by="SNP" )
intergen.all <- rbind(intergen.neut.merge, intergen.sel.merge)


# folded
ns.neut.merge$fold_bin <- cut(ns.neut.merge$folded_af, breaks=seq(0,0.5, 0.05))
ns.sel.merge$fold_bin <- cut(ns.sel.merge$folded_af, breaks=seq(0,0.5, 0.05))
ns.all$fold_bin <- cut(ns.all$folded_af, breaks=seq(0,0.5, 0.05))

syn.neut.merge$fold_bin <- cut(syn.neut.merge$folded_af, breaks=seq(0,0.5, 0.05))
syn.sel.merge$fold_bin <- cut(syn.sel.merge$folded_af, breaks=seq(0,0.5, 0.05))
syn.all$fold_bin <- cut(syn.all$folded_af, breaks=seq(0,0.5, 0.05))

intron.neut.merge$fold_bin <- cut(intron.neut.merge$folded_af, breaks=seq(0,0.5, 0.05))
intron.sel.merge$fold_bin <- cut(intron.sel.merge$folded_af, breaks=seq(0,0.5, 0.05))
intron.all$fold_bin <- cut(intron.all$folded_af, breaks=seq(0,0.5, 0.05))

intergen.neut.merge$fold_bin <- cut(intergen.neut.merge$folded_af, breaks=seq(0,0.5, 0.05))
intergen.sel.merge$fold_bin <- cut(intergen.sel.merge$folded_af, breaks=seq(0,0.5, 0.05))
intergen.all$fold_bin <- cut(intergen.all$folded_af, breaks=seq(0,0.5, 0.05))

# unfolded
ns.neut.merge$unfold_bin <- cut(ns.neut.merge$af_out, breaks=seq(0,1, 0.05))
ns.sel.merge$unfold_bin <- cut(ns.sel.merge$af_out, breaks=seq(0,1, 0.05))
ns.all$unfold_bin <- cut(ns.all$af_out, breaks=seq(0,1, 0.05))

syn.neut.merge$unfold_bin <- cut(syn.neut.merge$af_out, breaks=seq(0,1, 0.05))
syn.sel.merge$unfold_bin <- cut(syn.sel.merge$af_out, breaks=seq(0,1, 0.05))
syn.all$unfold_bin <- cut(syn.all$af_out, breaks=seq(0,1, 0.05))

intron.neut.merge$unfold_bin <- cut(intron.neut.merge$af_out, breaks=seq(0,1, 0.05))
intron.sel.merge$unfold_bin <- cut(intron.sel.merge$af_out, breaks=seq(0,1, 0.05))
intron.all$unfold_bin <- cut(intron.all$af_out, breaks=seq(0,1, 0.05))

intergen.neut.merge$unfold_bin <- cut(intergen.neut.merge$af_out, breaks=seq(0,1, 0.05))
intergen.sel.merge$unfold_bin <- cut(intergen.sel.merge$af_out, breaks=seq(0,1, 0.05))
intergen.all$unfold_bin <- cut(intergen.all$af_out, breaks=seq(0,1, 0.05))

### folded #####
bin <- c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5")

# make barplot of neutral vs selected
ns_neut_foldfreq <- table(ns.neut.merge$fold_bin)/sum(table(ns.neut.merge$fold_bin))
ns_sel_foldfreq <- table(ns.sel.merge$fold_bin)/sum(table(ns.sel.merge$fold_bin))
ns_all_foldfreq <- table(ns.all$fold_bin)/sum(table(ns.all$fold_bin))

syn_neut_foldfreq <- table(syn.neut.merge$fold_bin)/sum(table(syn.neut.merge$fold_bin))
syn_sel_foldfreq <- table(syn.sel.merge$fold_bin)/sum(table(syn.sel.merge$fold_bin))
syn_all_foldfreq <- table(syn.all$fold_bin)/sum(table(syn.all$fold_bin))

intron_neut_foldfreq <- table(intron.neut.merge$fold_bin)/sum(table(intron.neut.merge$fold_bin))
intron_sel_foldfreq <- table(intron.sel.merge$fold_bin)/sum(table(intron.sel.merge$fold_bin))
intron_all_foldfreq <- table(intron.all$fold_bin)/sum(table(intron.all$fold_bin))

intergen_neut_foldfreq <- table(intergen.neut.merge$fold_bin)/sum(table(intergen.neut.merge$fold_bin))
intergen_sel_foldfreq <- table(intergen.sel.merge$fold_bin)/sum(table(intergen.sel.merge$fold_bin))
intergen_all_foldfreq <- table(intergen.all$fold_bin)/sum(table(intergen.all$fold_bin))


# merge into tables for plotting
ns.sel.dat <- data.frame(class=rep("selected", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_sel_foldfreq))
ns.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_neut_foldfreq))
ns.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_all_foldfreq))
ns.fq.dat <- rbind(ns.sel.dat, ns.neut.dat, ns.all.dat)

syn.sel.dat <- data.frame(class=rep("selected", length(bin)), bin=bin,
            allele_frequency=as.vector(syn_sel_foldfreq))
syn.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(syn_neut_foldfreq))
syn.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(syn_all_foldfreq))
syn.fq.dat <- rbind(syn.sel.dat, syn.neut.dat, syn.all.dat)

intron.sel.dat <- data.frame(class=rep("selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intron_sel_foldfreq))
intron.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(intron_neut_foldfreq))
intron.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(intron_all_foldfreq))
intron.fq.dat <- rbind(intron.sel.dat, intron.neut.dat, intron.all.dat)

intergen.sel.dat <- data.frame(class=rep("selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intergen_sel_foldfreq))
intergen.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(intergen_neut_foldfreq))
intergen.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(intergen_all_foldfreq))
intergen.fq.dat <- rbind(intergen.sel.dat, intergen.neut.dat, intergen.all.dat)


a <- ggplot(data=ns.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("non-synonymous")+
  theme(legend.title=element_blank())

b <- ggplot(data=syn.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("synonymous")+
  theme(legend.title=element_blank())

c <- ggplot(data=intron.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("intron")+
  theme(legend.title=element_blank())

d <- ggplot(data=intergen.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("intergenic")+
  theme(legend.title=element_blank())


pdf("~/urchin_af/figures/AF_categories_folded_af.pdf", height=5, width=10)

grid.arrange(a, b, c, d, ncol=2)

dev.off()

####### Unfolded #########

bin <- c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25",
    "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5",
    "0.5-0.55", "0.55-0.6", "0.6-0.65", "0.65-0.7", "0.7-0.75",
    "0.75-0.8", "0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-1")

# make barplot of neutral vs selected
ns_neut_unfoldfreq <- table(ns.neut.merge$unfold_bin)/sum(table(ns.neut.merge$unfold_bin))
ns_sel_unfoldfreq <- table(ns.sel.merge$unfold_bin)/sum(table(ns.sel.merge$unfold_bin))
ns_all_unfoldfreq <- table(ns.all$unfold_bin)/sum(table(ns.all$unfold_bin))

syn_neut_unfoldfreq <- table(syn.neut.merge$unfold_bin)/sum(table(syn.neut.merge$unfold_bin))
syn_sel_unfoldfreq <- table(syn.sel.merge$unfold_bin)/sum(table(syn.sel.merge$unfold_bin))
syn_all_unfoldfreq <- table(syn.all$unfold_bin)/sum(table(syn.all$unfold_bin))

intron_neut_unfoldfreq <- table(intron.neut.merge$unfold_bin)/sum(table(intron.neut.merge$unfold_bin))
intron_sel_unfoldfreq <- table(intron.sel.merge$unfold_bin)/sum(table(intron.sel.merge$unfold_bin))
intron_all_unfoldfreq <- table(intron.all$unfold_bin)/sum(table(intron.all$unfold_bin))

intergen_neut_unfoldfreq <- table(intergen.neut.merge$unfold_bin)/sum(table(intergen.neut.merge$unfold_bin))
intergen_sel_unfoldfreq <- table(intergen.sel.merge$unfold_bin)/sum(table(intergen.sel.merge$unfold_bin))
intergen_all_unfoldfreq <- table(intergen.all$unfold_bin)/sum(table(intergen.all$unfold_bin))

# merge into tables for plotting
ns.sel.dat <- data.frame(class=rep("selected", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_sel_unfoldfreq))
ns.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_neut_unfoldfreq))
ns.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(ns_all_unfoldfreq))
ns.fq.dat <- rbind(ns.sel.dat, ns.neut.dat, ns.all.dat)

syn.sel.dat <- data.frame(class=rep("selected", length(bin)), bin=bin,
            allele_frequency=as.vector(syn_sel_unfoldfreq))
syn.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(syn_neut_unfoldfreq))
syn.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(syn_all_unfoldfreq))
syn.fq.dat <- rbind(syn.sel.dat, syn.neut.dat, syn.all.dat)

intron.sel.dat <- data.frame(class=rep("selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intron_sel_unfoldfreq))
intron.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(intron_neut_unfoldfreq))
intron.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(intron_all_unfoldfreq))
intron.fq.dat <- rbind(intron.sel.dat, intron.neut.dat, intron.all.dat)

intergen.sel.dat <- data.frame(class=rep("selected", length(bin)), bin=bin,
            allele_frequency=as.vector(intergen_sel_unfoldfreq))
intergen.neut.dat <- data.frame(class=rep("neutral", length(bin)),bin=bin,
            allele_frequency=as.vector(intergen_neut_unfoldfreq))
intergen.all.dat <- data.frame(class=rep("all variants", length(bin)),bin=bin,
            allele_frequency=as.vector(intergen_all_unfoldfreq))
intergen.fq.dat <- rbind(intergen.sel.dat, intergen.neut.dat, intergen.all.dat)



a <- ggplot(data=ns.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("non-synonymous")+
  theme(legend.title=element_blank())

b <- ggplot(data=syn.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("synonymous")+
  theme(legend.title=element_blank())

c <- ggplot(data=intron.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("intron")+
  theme(legend.title=element_blank())

d <- ggplot(data=intergen.fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("intergenic")+
  theme(legend.title=element_blank())

png("~/urchin_af/figures/AF_categories_unfolded_bin.png", height=7, res=300, units="in", width=10)
grid.arrange(a, b, c, d, ncol=2)
dev.off()


####### density plots

selected <- mydata[which(mydata$pH_sig==TRUE),]
neutral <- mydata[which(mydata$pH_sig==FALSE),]

# NS density plot

png("~/urchin_af/figures/AF_categories_unfolded_af.png", height=7, width=10, res=300, units="in")

par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 2.1))

plot(density(ns.neut.merge$af_out), ylim=c(0, 2.7), xlim=c(-0.05,1.05), lwd=3,
    xlab="allele frequency", col="blue", main="")
lines(density(ns.sel.merge$af_out), col="red", lwd=3)
#lines(density(ns.all$af_out), col="black", lwd=3)
lines(density(selected$af_out), col="red", lty=2, lwd=3)
#lines(density(neutral$af_out), col="black", lty=2, lwd=3)
title(main="non-synonymous")
#legend("topright", c("NS-neutral", "NS-selected", "NS-all", "Genome-wide: selected", "Genome-wide: neutral"),
 #   lty=c(1, 1, 1, 2, 2), col=c("blue", "red", "black", "red", "black"), lwd=2 )

plot(density(syn.neut.merge$af_out), ylim=c(0, 2.7), xlim=c(-0.05,1.05), lwd=3,
    xlab="allele frequency", col="blue", main="")
lines(density(syn.sel.merge$af_out), col="red", lwd=3)
#lines(density(syn.all$af_out), col="black", lwd=3)
lines(density(selected$af_out), col="red", lty=2, lwd=3)
#lines(density(neutral$af_out), col="black", lty=2, lwd=3)
title(main="synonymous")

plot(density(intron.neut.merge$af_out), ylim=c(0, 2.7), xlim=c(-0.05,1.05), lwd=3,
    xlab="allele frequency", col="blue", main="")
lines(density(intron.sel.merge$af_out), col="red", lwd=3)
#lines(density(intron.all$af_out), col="black", lwd=3)
lines(density(selected$af_out), col="red", lty=2, lwd=3)
#lines(density(neutral$af_out), col="black", lty=2, lwd=3)
title(main="intron")

plot(density(intergen.neut.merge$af_out), ylim=c(0, 2.7), xlim=c(-0.05,1.05), lwd=3,
    xlab="allele frequency", col="blue", main="")
lines(density(intergen.sel.merge$af_out), col="red", lwd=3)
#lines(density(intergen.all$af_out), col="black", lwd=3)
lines(density(selected$af_out), col="red", lty=2, lwd=3)
#lines(density(neutral$af_out), col="black", lty=2, lwd=3)
title(main="intergenic")


# add legend. basically overlaying new empty plot
par(fig = c(0.025, 1, 0, 1), oma = c(0, 0, .15, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0,0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("center", c("neutral", "selected", "Genome-wide: selected"),
    lty=c(1, 1, 2), col=c("blue", "red", "red"), lwd=2 )

dev.off()
