
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
dat$sig[which(dat$control_qval < 0.01)] <- TRUE
dat$sig[which(dat$PVAL < 0.01)] <- FALSE

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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.q01_ctr.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.q01_ctr.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
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


write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF_genes.q01_ctr.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF.q01_ctr.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC;

myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,
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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC_genes.q01_ctr.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC.q01_ctr.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


