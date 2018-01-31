

# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.


# Edit these to match your data file names:
#input="sig_binary.table" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
#goAnnotations="GO.annotation" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
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
dat <- read.table("~/urchin_af/analysis/go_enrichment/cmh.masterGO.out",header=TRUE)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)

geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.masterGO.annotation")  

# set gene background
geneUniverse <- names(geneID2GO) 

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
    topNodes = 29)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
#sigGenes(myGOdata)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP_genes.out", 
    sig_out, 
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_BP.out", allRes, 
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF; 

myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
    topNodes = 16)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
    ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF_genes.out", 
    sig_out, 
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_MF.out", allRes, 
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC; 

myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  
    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
    topNodes = 8)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
    ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC_genes.out", 
    sig_out, 
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_CC.out", allRes, 
    col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


###################################
##### coding
###################################

# read in gene list
dat <- read.table("~/urchin_af/analysis/go_enrichment/cmh.codingGO.out",header=TRUE)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)

geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.codingGO.annotation")  

# set gene background
geneUniverse <- names(geneID2GO) 

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP)

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
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms)){
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_coding_BP_genes.out", 
  sig_out, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_coding_BP.out", allRes, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF; 
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
  topNodes = 7)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_coding_MF_genes.out", 
  sig_out, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_coding_MF.out", allRes, 
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
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_coding_CC_genes.out", 
  sig_out, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_coding_CC.out", allRes, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")


###################################
##### NON-coding
###################################

# read in gene list
dat <- read.table("~/urchin_af/analysis/go_enrichment/cmh.NONcodingGO.out",header=TRUE)
i <- sapply(dat, is.factor)
dat[i] <- lapply(dat[i], as.character)

geneID2GO <- readMappings(file = "~/urchin_af/analysis/go_enrichment/topGO.NONcodingGO.annotation")  

# set gene background
geneUniverse <- names(geneID2GO) 

genesOfInterest <- dat[which(dat$sig == TRUE),]
genesOfInterest <- as.character(genesOfInterest$SNP)

#show genes of interest in universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
  topNodes = 55)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
sigGenes(myGOdata)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_NONcoding_BP_genes.out", 
  sig_out, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_NONcoding_BP.out", allRes, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for molecular function: MF; 
myGOdata <- new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
  topNodes = 42)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_NONcoding_MF_genes.out", 
  sig_out, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_NONcoding_MF.out", allRes, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

# rerun for Cellular Components: CC; 
myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  
  annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultWeight <- runTest(myGOdata, algorithm="weight", statistic="fisher") 
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher") 

allRes <- GenTable(myGOdata, classicFisher = resultClassic,  weight = resultWeight, orderBy = "weight", ranksOf = "weight", 
  topNodes = 12)

# find genes associated with sig go terms
myterms = allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
go_out <- as.data.frame(matrix(nrow=length(myterms), ncol=4))
sig_out <- as.data.frame(matrix(nrow=0, 
  ncol=(ncol(dat)+1)))
colnames(go_out) <- c("GO_term","gene_name", "WHL", "SPU")
for (i in 1:length(myterms))
   {
       mygenesforterm <- myterms[i]
       mygenesforterm <- mygenes[myterms][[i]]
       #make empty strings to append to
       whl_out <- c()
       spu_out <- c()
       gene_names <- c()

       for(z in 1:length(mygenesforterm)){
          out.1 <- dat[which(dat$SNP == mygenesforterm[z]),]
          out.2 <- cbind(myterms[i], out.1)
          sig_out <- rbind(sig_out, out.2) 
       }

     }

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_NONcoding_CC_genes.out", 
  sig_out, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")

write.table(file = "~/urchin_af/analysis/go_enrichment/go_enrichment_NONcoding_CC.out", allRes, 
  col.names=TRUE, row.names=FALSE, quote=FALSE,sep="\t")



