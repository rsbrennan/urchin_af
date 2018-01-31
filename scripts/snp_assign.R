# snp assignments

library(stringr)
library(ggplot2)
library(gridExtra)
library(MASS)


dat <- read.table("~/urchin_af/variants/urchin_ann.vcf", stringsAsFactors=FALSE)

cmh <- read.table("~/urchin_af/analysis/cmh.master.out", header=TRUE)

# it is possible that snps receive multiple annotations
# I think in this case, it is best to take the most "severe"
# order is: HIGH MODERATE LOW MODIFIER

# loop over entire data set. if more than one annotation
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

# make empty dataframe
out <- dat

for (i in 1:length(eff)){
    if(length(eff[[i]]) == 1){
        out$V8[i] <- eff[[i]]
    }
    else{
        # pull out the change of each variant
        change <- lapply(strsplit(eff[[i]], split="|", fixed=TRUE), '[[', 3)
        #check if they all match. if so, just take the first one
        if(length(unique(change)) == 1){
            out$V8[i] <- eff[[i]][1]
        }
        #if they don't all match, take the one with highest change
        else{
            change <- gsub("MODIFIER", "1", change)
            change <- gsub("LOW", "2", change)
            change <- gsub("MODERATE", "3", change)
            change <- gsub("HIGH", "4", change)
            ef_max <- (max(change))
            pick_max <- which(change == ef_max)
            if(length(pick_max) == 1){
                out$V8[i] <- eff[[i]][pick_max]
            }
            else{
                out$V8[i] <- eff[[i]][pick_max[1]]
            }
        }
    }
    if(i%%1000 == 0){print(i)}
}

out$SNP <- paste(out$V1, out$V2, sep=":")

out_ann <- data.frame(SNP=out$SNP, ANN=out$V8)
# split annotation
ann_split <- (str_split_fixed(out$V8, "\\|", n=16))
out_ann <-(cbind(out$SNP, ann_split))

colnames(out_ann) <- c("SNP","Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID", "Feature_Type","Feature_ID", "Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos",  "cDNA.length","CDS.pos","Distance","WARNINGS")

new <- merge(cmh,out_ann, by="SNP")

write.table(file="~/urchin_af/analysis/cmh.annotations.out", new, 
    col.name=TRUE, quote=FALSE, row.name=FALSE )


#### Repeat for custom annotation using Transcriptome
dat <- read.table("~/urchin_af/variants/urchin_ann_custom.vcf", stringsAsFactors=FALSE)

cmh <- read.table("~/urchin_af/analysis/cmh.master.out", header=TRUE)

# it is possible that snps receive multiple annotations
# I think in this case, it is best to take the most "severe"
# order is: HIGH MODERATE LOW MODIFIER

# loop over entire data set. if more than one annotation
eff <- strsplit(as.character(dat$V8), split=",", fixed=TRUE)

# make  dataframe
out <- dat

for (i in 1:length(eff)){
    if(length(eff[[i]]) == 1){
        out$V8[i] <- eff[[i]]
    }
    else{
        # pull out the change of each variant
        change <- lapply(strsplit(eff[[i]], split="|", fixed=TRUE), '[[', 3)
        #check if they all match. if so, just take the first one
        if(length(unique(change)) == 1){
            out$V8[i] <- eff[[i]][1]
        }
        #if they don't all match, take the one with highest change
        else{
            change <- gsub("MODIFIER", "1", change)
            change <- gsub("LOW", "2", change)
            change <- gsub("MODERATE", "3", change)
            change <- gsub("HIGH", "4", change)
            ef_max <- (max(change))
            pick_max <- which(change == ef_max)
            if(length(pick_max) == 1){
                out$V8[i] <- eff[[i]][pick_max]
            }
            else{
                out$V8[i] <- eff[[i]][pick_max[1]]
            }
        }
    }
    if(i%%5000 == 0){print(i)}
}

out$SNP <- paste(out$V1, out$V2, sep=":")

out_ann <- data.frame(SNP=out$SNP, ANN=out$V8)
# split annotation
ann_split <- (str_split_fixed(out$V8, "\\|", n=16))
out_ann <-(cbind(out$SNP, ann_split))

colnames(out_ann) <- c("SNP","Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID", "Feature_Type","Feature_ID", "Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos",  "cDNA.length","CDS.pos","Distance","WARNINGS")

new_cust <- merge(cmh,out_ann, by="SNP")

write.table(file="~/urchin_af/analysis/cmh.annotations_custom.out", new_cust, 
    col.name=TRUE, quote=FALSE, row.name=FALSE )

####
# assign functional categories
####

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

ens <- as.data.frame(matrix(ncol=3, nrow=5))
colnames(ens) <- c("gff", "SNP_class", "count")
ens$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous", "3_prime_UTR")
ens$gff <- "ensembl"
ens$count[1] <- length(which(new$Annotation == "downstream_gene_variant" | 
    new$Annotation == "intergenic_region" |
    new$Annotation == "upstream_gene_variant"))
ens$count[2] <- length(which(new$Annotation == "intron_variant" | 
    new$Annotation == "splice_region_variant&intron_variant"| 
    new$Annotation == "splice_acceptor_variant&intron_variant"| 
    new$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"| 
    new$Annotation == "splice_region_variant"| 
    new$Annotation == "splice_region_variant&intron_variant" |
    new$Annotation == "splice_region_variant&synonymous_variant"  ))
ens$count[3] <- length(which(new$Annotation == "synonymous_variant" | 
    new$Annotation == "stop_retained_variant"))
ens$count[4] <- length(which(new$Annotation == "missense_variant" | 
    new$Annotation == "stop_gained" | 
    new$Annotation == "stop_lost"| 
    new$Annotation == "missense_variant&splice_region_variant"| 
    new$Annotation == "start_lost"| 
    new$Annotation == "stop_gained&splice_region_variant"| 
    new$Annotation == "stop_lost&splice_region_variant"| 
    new$Annotation == "initiator_codon_variant"))
ens$count[5] <- length(which(new$Annotation == "3_prime_UTR_variant"))
ens$proportion <- ens$count/sum(ens$count)


# rerun for transcritpome gff
cust <- as.data.frame(matrix(ncol=3, nrow=5))
colnames(cust) <- c("gff", "SNP_class", "count")
cust$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous", "3_prime_UTR")
cust$gff <- "Transcriptome"
cust$count[1] <- length(which(new_cust$Annotation == "downstream_gene_variant" | 
    new_cust$Annotation == "intergenic_region" |
    new_cust$Annotation == "upstream_gene_variant"))
cust$count[2] <- length(which(new_cust$Annotation == "intron_variant" | 
    new_cust$Annotation == "splice_region_variant&intron_variant"| 
    new_cust$Annotation == "splice_acceptor_variant&intron_variant"| 
    new_cust$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"| 
    new_cust$Annotation == "splice_region_variant"| 
    new_cust$Annotation == "splice_region_variant&intron_variant" |
    new_cust$Annotation == "splice_region_variant&synonymous_variant"  ))
cust$count[3] <- length(which(new_cust$Annotation == "synonymous_variant" | 
    new_cust$Annotation == "stop_retained_variant"))
cust$count[4] <- length(which(new$Annotation == "missense_variant" | 
    new_cust$Annotation == "stop_gained" | 
    new_cust$Annotation == "stop_lost"| 
    new_cust$Annotation == "missense_variant&splice_region_variant"| 
    new_cust$Annotation == "start_lost"| 
    new_cust$Annotation == "stop_gained&splice_region_variant"| 
    new_cust$Annotation == "stop_lost&splice_region_variant"| 
    new_cust$Annotation == "initiator_codon_variant"))
cust$count[5] <- length(which(new_cust$Annotation == "3_prime_UTR_variant"))
cust$proportion <- cust$count/sum(cust$count)

#bind the two together
dat_count <- rbind(ens, cust)

pdf("~/urchin_af/figures/snp_class.pdf", height=7, width=7)

a <- ggplot(data=dat_count, aes(x=SNP_class, y=proportion, fill=gff)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3'))
a

dev.off()

############# 
# split into selected and non
############# 

ens.sel <- as.data.frame(matrix(ncol=3, nrow=5))
colnames(ens.sel) <- c("gff", "SNP_class", "count")
ens.sel$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous", "3_prime_UTR")
ens.sel$gff <- "selected"
new.sel <- new[which(new$sig ==TRUE),]

ens.sel$count[1] <- length(which(new.sel$Annotation == "downstream_gene_variant" | 
    new.sel$Annotation == "intergenic_region" |
    new.sel$Annotation == "upstream_gene_variant"))
ens.sel$count[2] <- length(which(new.sel$Annotation == "intron_variant" | 
    new.sel$Annotation == "splice_region_variant&intron_variant"| 
    new.sel$Annotation == "splice_acceptor_variant&intron_variant"| 
    new.sel$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"| 
    new.sel$Annotation == "splice_region_variant"| 
    new.sel$Annotation == "splice_region_variant&intron_variant" |
    new.sel$Annotation == "splice_region_variant&synonymous_variant"  ))
ens.sel$count[3] <- length(which(new.sel$Annotation == "synonymous_variant" | 
    new.sel$Annotation == "stop_retained_variant"))
ens.sel$count[4] <- length(which(new.sel$Annotation == "missense_variant" | 
    new.sel$Annotation == "stop_gained" | 
    new.sel$Annotation == "stop_lost"| 
    new.sel$Annotation == "missense_variant&splice_region_variant"| 
    new.sel$Annotation == "start_lost"| 
    new.sel$Annotation == "stop_gained&splice_region_variant"| 
    new.sel$Annotation == "stop_lost&splice_region_variant"| 
    new.sel$Annotation == "initiator_codon_variant"))
ens.sel$count[5] <- length(which(new.sel$Annotation == "3_prime_UTR_variant"))
ens.sel$proportion <- ens.sel$count/sum(ens.sel$count)

#neutral loci
ens.neut <- as.data.frame(matrix(ncol=3, nrow=5))
colnames(ens.neut) <- c("gff", "SNP_class", "count")
ens.neut$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous", "3_prime_UTR")
ens.neut$gff <- "neutral"
new.neut <- new[which(new$sig ==FALSE),]

ens.neut$count[1] <- length(which(new.neut$Annotation == "downstream_gene_variant" | 
    new.neut$Annotation == "intergenic_region" |
    new.neut$Annotation == "upstream_gene_variant"))
ens.neut$count[2] <- length(which(new.neut$Annotation == "intron_variant" | 
    new.neut$Annotation == "splice_region_variant&intron_variant"| 
    new.neut$Annotation == "splice_acceptor_variant&intron_variant"| 
    new.neut$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"| 
    new.neut$Annotation == "splice_region_variant"| 
    new.neut$Annotation == "splice_region_variant&intron_variant" |
    new.neut$Annotation == "splice_region_variant&synonymous_variant"  ))
ens.neut$count[3] <- length(which(new.neut$Annotation == "synonymous_variant" | 
    new.neut$Annotation == "stop_retained_variant"))
ens.neut$count[4] <- length(which(new.neut$Annotation == "missense_variant" | 
    new.neut$Annotation == "stop_gained" | 
    new.neut$Annotation == "stop_lost"| 
    new.neut$Annotation == "missense_variant&splice_region_variant"| 
    new.neut$Annotation == "start_lost"| 
    new.neut$Annotation == "stop_gained&splice_region_variant"| 
    new.neut$Annotation == "stop_lost&splice_region_variant"| 
    new.neut$Annotation == "initiator_codon_variant"))
ens.neut$count[5] <- length(which(new.neut$Annotation == "3_prime_UTR_variant"))
ens.neut$proportion <- ens.neut$count/sum(ens.neut$count)

ens_all <- rbind(ens.neut, ens.sel)

a <- ggplot(data=ens_all, aes(x=SNP_class, y=proportion, fill=gff)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3')) +
  ggtitle("Ensembl gff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

a


# custom gff
cust.sel <- as.data.frame(matrix(ncol=3, nrow=5))
colnames(cust.sel) <- c("gff", "SNP_class", "count")
cust.sel$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous", "3_prime_UTR")
cust.sel$gff <- "selected"
new.sel <- new_cust[which(new_cust$sig ==TRUE),]

cust.sel$count[1] <- length(which(new.sel$Annotation == "downstream_gene_variant" | 
    new.sel$Annotation == "intergenic_region" |
    new.sel$Annotation == "upstream_gene_variant"))
cust.sel$count[2] <- length(which(new.sel$Annotation == "intron_variant" | 
    new.sel$Annotation == "splice_region_variant&intron_variant"| 
    new.sel$Annotation == "splice_acceptor_variant&intron_variant"| 
    new.sel$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"| 
    new.sel$Annotation == "splice_region_variant"| 
    new.sel$Annotation == "splice_region_variant&intron_variant" |
    new.sel$Annotation == "splice_region_variant&synonymous_variant"  ))
cust.sel$count[3] <- length(which(new.sel$Annotation == "synonymous_variant" | 
    new.sel$Annotation == "stop_retained_variant"))
cust.sel$count[4] <- length(which(new.sel$Annotation == "missense_variant" | 
    new.sel$Annotation == "stop_gained" | 
    new.sel$Annotation == "stop_lost"| 
    new.sel$Annotation == "missense_variant&splice_region_variant"| 
    new.sel$Annotation == "start_lost"| 
    new.sel$Annotation == "stop_gained&splice_region_variant"| 
    new.sel$Annotation == "stop_lost&splice_region_variant"| 
    new.sel$Annotation == "initiator_codon_variant"))
cust.sel$count[5] <- length(which(new.sel$Annotation == "3_prime_UTR_variant"))
cust.sel$proportion <- cust.sel$count/sum(cust.sel$count)

cust.neut <- as.data.frame(matrix(ncol=3, nrow=5))
colnames(cust.neut) <- c("gff", "SNP_class", "count")
cust.neut$SNP_class <- c("intergenic", "intron", "synonymous", "non-synonymous", "3_prime_UTR")
cust.neut$gff <- "neutral"
new.neut <- new_cust[which(new_cust$sig ==FALSE),]

cust.neut$count[1] <- length(which(new.neut$Annotation == "downstream_gene_variant" | 
    new.neut$Annotation == "intergenic_region" |
    new.neut$Annotation == "upstream_gene_variant"))
cust.neut$count[2] <- length(which(new.neut$Annotation == "intron_variant" | 
    new.neut$Annotation == "splice_region_variant&intron_variant"| 
    new.neut$Annotation == "splice_acceptor_variant&intron_variant"| 
    new.neut$Annotation == "splice_donor_variant&splice_region_variant&intron_variant"| 
    new.neut$Annotation == "splice_region_variant"| 
    new.neut$Annotation == "splice_region_variant&intron_variant" |
    new.neut$Annotation == "splice_region_variant&synonymous_variant"  ))
cust.neut$count[3] <- length(which(new.neut$Annotation == "synonymous_variant" | 
    new.neut$Annotation == "stop_retained_variant"))
cust.neut$count[4] <- length(which(new.neut$Annotation == "missense_variant" | 
    new.neut$Annotation == "stop_gained" | 
    new.neut$Annotation == "stop_lost"| 
    new.neut$Annotation == "missense_variant&splice_region_variant"| 
    new.neut$Annotation == "start_lost"| 
    new.neut$Annotation == "stop_gained&splice_region_variant"| 
    new.neut$Annotation == "stop_lost&splice_region_variant"| 
    new.neut$Annotation == "initiator_codon_variant"))
cust.neut$count[5] <- length(which(new.neut$Annotation == "3_prime_UTR_variant"))
cust.neut$proportion <- cust.neut$count/sum(cust.neut$count)

cust_all <- rbind(cust.neut, cust.sel)

b <- ggplot(data=cust_all, aes(x=SNP_class, y=proportion, fill=gff)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3'))+
  ggtitle("Transcriptome gff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

b

pdf("~/urchin_af/figures/snp_eff.pdf", height=5, width=10)
grid.arrange(a, b, ncol=2)
dev.off()

# make contingency table

r1 <- c(24280, 40542-24280)
r2 <- c(506, 780-506)

din <- matrix(c(r1, r2),
                nrow=2,
                byrow=TRUE)
rownames(din) <- c("neutral", "selected")
colnames(din) <- c("non-synonymous", "all") 

chisq.test(din,correct=FALSE)    

r1 <- c(4241, 40501-4241)
r2 <- c(94, 780-94)

din <- matrix(c(r1, r2),
                nrow=2,
                byrow=TRUE)
rownames(din) <- c("neutral", "selected")
colnames(din) <- c("non-synonymous", "all") 

chisq.test(din,correct=FALSE)   

########
### compare coding vs non
########

r1 <- c((14776+7647), (12439+4241+1398))
r2 <- c((281+142),(238+94+25))

din <- matrix(c(r1, r2),
                nrow=2,
                byrow=TRUE)
rownames(din) <- c("neutral", "selected")
colnames(din) <- c("coding", "non-coding")                      
chisq.test(din,correct=FALSE)    

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

## ensembl gene models

new.sel <- new[which(new$sig ==TRUE),]
ens.sel <-  new.sel[which(new.sel$Annotation == "missense_variant" | 
    new.sel$Annotation == "stop_gained" | 
    new.sel$Annotation == "stop_lost"| 
    new.sel$Annotation == "missense_variant&splice_region_variant"| 
    new.sel$Annotation == "start_lost"| 
    new.sel$Annotation == "stop_gained&splice_region_variant"| 
    new.sel$Annotation == "stop_lost&splice_region_variant"| 
    new.sel$Annotation == "initiator_codon_variant"),]


#neutral loci
new.neut <- new[which(new$sig ==FALSE),]
ens.neut <- new.neut[which(new.neut$Annotation == "missense_variant" | 
    new.neut$Annotation == "stop_gained" | 
    new.neut$Annotation == "stop_lost"| 
    new.neut$Annotation == "missense_variant&splice_region_variant"| 
    new.neut$Annotation == "start_lost"| 
    new.neut$Annotation == "stop_gained&splice_region_variant"| 
    new.neut$Annotation == "stop_lost&splice_region_variant"| 
    new.neut$Annotation == "initiator_codon_variant"),]

# transcriptome gene models

new.sel <- new_cust[which(new_cust$sig ==TRUE),]
cust.sel <- new.sel[which(new.sel$Annotation == "missense_variant" | 
    new.sel$Annotation == "stop_gained" | 
    new.sel$Annotation == "stop_lost"| 
    new.sel$Annotation == "missense_variant&splice_region_variant"| 
    new.sel$Annotation == "start_lost"| 
    new.sel$Annotation == "stop_gained&splice_region_variant"| 
    new.sel$Annotation == "stop_lost&splice_region_variant"| 
    new.sel$Annotation == "initiator_codon_variant"),]

new.neut <- new_cust[which(new_cust$sig ==FALSE),]

cust.neut <- new.neut[which(new.neut$Annotation == "missense_variant" | 
    new.neut$Annotation == "stop_gained" | 
    new.neut$Annotation == "stop_lost"| 
    new.neut$Annotation == "missense_variant&splice_region_variant"| 
    new.neut$Annotation == "start_lost"| 
    new.neut$Annotation == "stop_gained&splice_region_variant"| 
    new.neut$Annotation == "stop_lost&splice_region_variant"| 
    new.neut$Annotation == "initiator_codon_variant"),]

# read in allele frequency data

mydata <- read.table("~/urchin_af/analysis/adaptive_allelefreq.txt", header=TRUE)
mydata$SNP <- paste(mydata$CHROM, mydata$POS, sep=":")

# fold af

mydata$folded_af = unlist(lapply(mydata$D1_8_mean,function(x)  
          ifelse(x > 0.5, (1-x), x)), use.names=FALSE)
# unfolded af
# these allele freqs are col: af_out

# merge datasets:

cust.neut.merge <- merge(mydata,cust.neut, by="SNP" )
cust.sel.merge <- merge(mydata,cust.sel, by="SNP" )
ens.neut.merge <- merge(mydata,ens.neut, by="SNP" )
ens.sel.merge <- merge(mydata,ens.sel, by="SNP" )

cust.all <- rbind(cust.neut.merge,cust.sel.merge)
ens.all <- rbind(ens.neut.merge,ens.sel.merge)

# folded

cust.neut.merge$fold_bin <- cut(cust.neut.merge$folded_af, breaks=seq(0,0.5, 0.05))
cust.sel.merge$fold_bin <- cut(cust.sel.merge$folded_af, breaks=seq(0,0.5, 0.05))
ens.neut.merge$fold_bin <- cut(ens.neut.merge$folded_af, breaks=seq(0,0.5, 0.05))
ens.sel.merge$fold_bin <- cut(ens.sel.merge$folded_af, breaks=seq(0,0.5, 0.05))

ens.all$fold_bin <- cut(ens.all$folded_af, breaks=seq(0,0.5, 0.05))
cust.all$fold_bin <- cut(cust.all$folded_af, breaks=seq(0,0.5, 0.05))

# unfolded

cust.neut.merge$unfold_bin <- cut(cust.neut.merge$af_out, breaks=seq(0,1, 0.05))
cust.sel.merge$unfold_bin <- cut(cust.sel.merge$af_out, breaks=seq(0,1, 0.05))
ens.neut.merge$unfold_bin <- cut(ens.neut.merge$af_out, breaks=seq(0,1, 0.05))
ens.sel.merge$unfold_bin <- cut(ens.sel.merge$af_out, breaks=seq(0,1, 0.05))

ens.all$unfold_bin <- cut(ens.all$af_out, breaks=seq(0,1, 0.05))
cust.all$unfold_bin <- cut(cust.all$af_out, breaks=seq(0,1, 0.05))

# make barplot of neutral vs selected transcriptom

bin <- c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5")

cust_neut_foldfreq <- table(cust.neut.merge$fold_bin)/sum(table(cust.neut.merge$fold_bin))
cust_sel_foldfreq <- table(cust.sel.merge$fold_bin)/sum(table(cust.sel.merge$fold_bin))
cust_all_foldfreq <- table(cust.all$fold_bin)/sum(table(cust.all$fold_bin))

sel.dat <- data.frame(class=rep("selected", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(cust_sel_foldfreq))
neut.dat <- data.frame(class=rep("neutral", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(cust_neut_foldfreq))
all.dat <- data.frame(class=rep("all variants", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(cust_all_foldfreq))


fq.dat <- rbind(sel.dat,neut.dat, all.dat)


b <- ggplot(data=fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  ggtitle("Transcriptome gff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b

# make barplot of neutral vs selected ensembl
ens_neut_foldfreq <- table(ens.neut.merge$fold_bin)/sum(table(ens.neut.merge$fold_bin))
ens_sel_foldfreq <- table(ens.sel.merge$fold_bin)/sum(table(ens.sel.merge$fold_bin))
ens_all_foldfreq <- table(ens.all$fold_bin)/sum(table(ens.all$fold_bin))

sel.dat <- data.frame(class=rep("selected", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(ens_sel_foldfreq))
neut.dat <- data.frame(class=rep("neutral", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(ens_neut_foldfreq))
all.dat <- data.frame(class=rep("all variants", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(ens_all_foldfreq))
fq.dat <- rbind(sel.dat,neut.dat, all.dat)


a <- ggplot(data=fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  ggtitle("Ensembl gff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#a

pdf("~/urchin_af/figures/NS_folded_af.pdf", height=5, width=10)
grid.arrange(a, b, ncol=2)
dev.off()

# rerun for unfolded allele frequency

# make barplot of neutral vs selected transcriptom

bin <- c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", 
    "0.25-0.3", "0.3-0.35", "0.35-0.4", "0.4-0.45", "0.45-0.5",
    "0.5-0.55", "0.55-0.6", "0.6-0.65", "0.65-0.7", "0.7-0.75", 
    "0.75-0.8", "0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-1")

cust_neut_unfoldfreq <- table(cust.neut.merge$unfold_bin)/sum(table(cust.neut.merge$unfold_bin))
cust_sel_unfoldfreq <- table(cust.sel.merge$unfold_bin)/sum(table(cust.sel.merge$unfold_bin))
cust_all_unfoldfreq <- table(cust.all$unfold_bin)/sum(table(cust.all$unfold_bin))

sel.dat <- data.frame(class=rep("selected", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(cust_sel_unfoldfreq))
neut.dat <- data.frame(class=rep("neutral", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(cust_neut_unfoldfreq))
all.dat <- data.frame(class=rep("all variants", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(cust_all_unfoldfreq))


fq.dat <- rbind(sel.dat,neut.dat, all.dat)


b <- ggplot(data=fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  ggtitle("Transcriptome gff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#b

# make barplot of neutral vs selected ensembl
ens_neut_unfoldfreq <- table(ens.neut.merge$unfold_bin)/sum(table(ens.neut.merge$unfold_bin))
ens_sel_unfoldfreq <- table(ens.sel.merge$unfold_bin)/sum(table(ens.sel.merge$unfold_bin))
ens_all_unfoldfreq <- table(ens.all$unfold_bin)/sum(table(ens.all$unfold_bin))

sel.dat <- data.frame(class=rep("selected", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(ens_sel_unfoldfreq))
neut.dat <- data.frame(class=rep("neutral", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(ens_neut_unfoldfreq))
all.dat <- data.frame(class=rep("all variants", length(bin)),
            bin=bin, 
            allele_frequency=as.vector(ens_all_unfoldfreq))
fq.dat <- rbind(sel.dat,neut.dat, all.dat)


a <- ggplot(data=fq.dat, aes(x=bin, y=allele_frequency, fill=class, group=class)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_manual(values=c('royalblue4','tomato3', 'grey38'))+
  ggtitle("Ensembl gff") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#a

pdf("~/urchin_af/figures/NS_unfolded_af.pdf", height=5, width=10)
grid.arrange(a, b, ncol=2)
dev.off()


####### density plots

selected <- mydata[which(mydata$pH_sig==TRUE),]
neutral <- mydata[which(mydata$pH_sig==FALSE),]

png("~/urchin_af/figures/NS_unfolded_af.png", height=5, width=10, res=300, units="in")
plot(density(ens.neut.merge$af_out), ylim=c(0, 2.7), xlim=c(-0.05,1.05), lwd=3, 
    xlab="allele frequency", col="blue", main="")
lines(density(ens.sel.merge$af_out), col="red", lwd=3)
lines(density(ens.all$af_out), col="black", lwd=3)
lines(density(selected$af_out), col="red", lty=2, lwd=3)
lines(density(neutral$af_out), col="black", lty=2, lwd=3)

legend("topright", c("NS-neutral", "NS-selected", "NS-all", "Genome-wide: selected", "Genome-wide: neutral"), 
    lty=c(1, 1, 1, 2, 2), col=c("blue", "red", "black", "red", "black"), lwd=3 )
dev.off()

