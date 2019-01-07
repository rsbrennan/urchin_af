
### GO enrichment

# infiles are: master genic exon intergenic synonymous non_syn intron non_coding

############################################################
#############  topGO
############################################################

library(topGO)

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.pH75_master_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig_pH75)
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.pH75.master.annotation")

#dat$sig <- FALSE
#dat$sig[which(dat$qval_pH75 < 0.001)] <- TRUE

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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.pH75.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.pH75.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

######################################################################
###################################
## pH 8.0
###################################
######################################################################

dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.pH80_master_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig_pH80)
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.pH80.master.annotation")

#dat$sig <- FALSE
#dat$sig[which(dat$qval_pH80 < 0.001)] <- TRUE
#dat$sig[which(dat$qval_pH75 < 0.001)] <- FALSE

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
    topNodes = 10)

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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.pH80.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.pH80.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


######################################################################
###################################
##  overlap
###################################
######################################################################

dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.overlap_master_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.overlap.master.annotation")

# set gene background
geneUniverse <- names(geneID2GO)

genesOfInterest <- dat[which(dat$overlap_sig == "TRUE"),]
genesOfInterest <- as.character(genesOfInterest$SNP_1)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher")
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight",
    topNodes = 19)

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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.overlap.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.overlap.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

