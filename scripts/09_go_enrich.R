
### GO enrichment

# infiles are: master genic exon intergenic synonymous non_syn intron non_coding

############################################################
#############  topGO
############################################################

library(topGO)

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.master_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.master.annotation")

dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == "TRUE"),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
    topNodes = 14)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.q01.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.q01.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
    topNodes = 14)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }


write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF_genes.q01.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF.q01.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;

myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
    topNodes = 4)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC_genes.q01.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC.q01.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


###################################
##### genic
###################################
# read in gene list

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.genic_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.genic.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 15)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_genic_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_genic_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 17)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_genic_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_genic_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 4)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_genic_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_genic_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


###################################
##### coding
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.exon_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.exon.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 11)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_exon_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_exon_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 11)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_exon_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_exon_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 7)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_exon_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_exon_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

###################################
##### intergenic
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.intergenic_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.intergenic.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intergenic_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intergenic_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intergenic_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intergenic_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intergenic_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intergenic_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

###################################
##### synonymous
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.synonymous_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.synonymous.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_synonymous_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_synonymous_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_synonymous_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_synonymous_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_synonymous_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_synonymous_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

###################################
##### non_syn
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.non_syn_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.non_syn.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_syn_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_syn_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_syn_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_syn_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_syn_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_syn_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

###################################
##### intron
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.intron_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.intron.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intron_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intron_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = 2)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intron_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intron_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intron_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_intron_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

###################################
##### non_coding
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.non_coding_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
dat$sig <- FALSE
dat$sig[which(dat$PVAL < 0.01)] <- TRUE
dat$sig[which(dat$control_qval < 0.01)] <- FALSE
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.non_coding.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_BP_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_BP.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_MF_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_MF.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
  topNodes = length(which(as.numeric(GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",topNodes = 10)$weight) < 0.05)))

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
sig_out <- as.data.frame(matrix(nrow=0,
  ncol=(ncol(dat)+1)))
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP_1 == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2)
       }
     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
