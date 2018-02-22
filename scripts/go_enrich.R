
### GO enrichment

# infiles are: master genic exon intergenic synonymous non_syn intron non_coding


# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.

# Edit these to match your data file names:
#input="pval.master.table" # two columns of comma-separated values: gene id, continuous measure of significance. 
  # To perform standard GO enrichment analysis based on 
  # Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
#goAnnotations="GO.master.annotation" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
#goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
#goDivision="MF" # either MF, or BP, or CC
#source("gomwu.functions.R")
#
#
## Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
#gomwuStats(input, goDatabase, goAnnotations, goDivision,
#   perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
#   largest=0.4,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
#   smallest=2,   # a GO category should contain at least this many genes to be considered
#   clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#)
## do not continue if the printout shows that no GO terms pass 10% FDR.
#
#
## Plotting results
#quartz()
#results=gomwuPlot(input,goAnnotations,goDivision,
##  absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
#   absValue=1,
#   level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
#   level2=0.05, # FDR cutoff to print in regular (not italic) font.
#   level3=0.01, # FDR cutoff to print in large bold font.
#   txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
#   treeHeight=0.5, # height of the hierarchical clustering tree
##  colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
#)
## manually rescale the plot so the tree matches the text
## if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.
#
## text representation of results, with actual adjusted p-values
#results

############################################################
#############  topGO
############################################################

library(topGO)

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.master_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.master.annotation")

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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF;

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,
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


write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF_genes.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF.txt", allRes,
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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC_genes.txt",
    sig_out,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC.txt", allRes,
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


###################################
##### genic
###################################

# read in gene list
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.genic_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
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
  topNodes = 12)

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
  topNodes = 8)

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
  topNodes = 3)

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
  topNodes =1)

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
  topNodes = 12)

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
  topNodes = 5)

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
dat <- read.csv("~/urchin_af/analysis/go_enrichment/cmh.synonymous_GO.out",header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat$sig <- as.logical(dat$sig)
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
  topNodes = 5)

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
  topNodes = 5)

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
  topNodes = 5)

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
  topNodes = 5)

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
  topNodes = 1)

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
  topNodes = 6)

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
  topNodes = 5)

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

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_CC_genes.txt",
  sig_out,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_non_coding_CC.txt", allRes,
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")
