###########################################################################
# Created 25 Feb 2015
# BioC 3.0

# DE analysis of Affymetrix Mouse Gene 2.0 ST arrays (pd.mogene.2.0.st)

# Update 9 Mar 2015 
# - find Ensebl IDs
# 11 Mar 2015
# - add Ensembl annotation 
# - filtering based on GC content
# 12 Mar 2015
# - distributions of expression stratified by GC content



###########################################################################

setwd("/home/Shared/data/array/Microarray_Edwin")

###########################################################################
#### create targets table with information about samples
###########################################################################

library(limma)


metadata <- read.table("metadata/Micro-array sample list.txt", sep = "\t")

targets <- data.frame(metadata , filename = list.files("CEL/", pattern="IA2" ,full.names = TRUE))

colnames(targets) <- c("Experiment", "SampleNr", "CellType", "FileName")

targets$ExperimentShort <- targets$Experiment
levels(targets$ExperimentShort) <- c("Bone marrow control"="control", "Kit control"="control", "Leukemia"="leukemia", "Leukemia after treatment" = "afterTreatment", "T cell control" = "control", "Thymocyte control" = "control")

targets$CellTypeShort <- targets$CellType
levels(targets$CellTypeShort) <- c("907" = "907", "907 - Post Dex" = "907", "B2M10" = "B2M10", "B2M2" = "B2M2", "B2M3" = "B2M3", "B2M3 Post dex" = "B2M3", "B3M3" = "B3M3", "B3M30" = "B3M30", "CD4 T cells spleen 1" = "CD4", "CD4 T cells spleen 2" = "CD4", "CD4 T cells spleen 3" = "CD4", "CD4+8+ DP Thymocytes 1" = "CD4+8+", "CD4+8+ DP Thymocytes 2" = "CD4+8+", "CD4+8+ DP Thymocytes 3" = "CD4+8+", "CD8 T cells spleen 1" = "CD8", "CD8 T cells spleen 2" = "CD8", "CD8 T cells spleen 3" = "CD8", "HeLa control" = "HeLa", "Whole bone marrow 1" = "wholeBoneMarrow", "Whole bone marrow 2" = "wholeBoneMarrow", "Whole bone marrow 3" = "wholeBoneMarrow")
  
targets$labels <- factor(paste(targets$ExperimentShort, targets$CellTypeShort, sep="_" ))

targets$groups <- targets$labels
levels(targets$groups)[grep(pattern = "leukemia", levels(targets$groups))] <- "leukemia"
levels(targets$groups)[grep(pattern = "afterTreatment", levels(targets$groups))] <- "afterTreatment"

nlevels(targets$groups)

targets$ctrlRep <- c(rep("", 12), rep(1:3, rep(4, 3)))


library(RColorBrewer)
# mypalette <- colorRampPalette(brewer.pal(11, "Spectral"))
# mypalette <- colorRampPalette(brewer.pal(9, "Set1"))
# mypalette <- function(n) {
#   hues = seq(15, 375, length=n+1)
#   hcl(h=hues, l=65, c=100)[1:n]
# }
# 
# mypalette <- colorRampPalette(c("brown", "red", "orange", "yellow", "green", "blue", "violet", "pink", "grey"))

# mypalette <- colorRampPalette(c("coral4", "firebrick2", "darkorange2", "goldenrod2", "forestgreen", "dodgerblue3", "darkviolet", "deeppink", "antiquewhite4"))

mypalette <- colorRampPalette(c("firebrick2", "darkorange2", "goldenrod2", "forestgreen", "dodgerblue3","darkviolet", "deeppink", "lavenderblush4"))


# p <- 12
# plot(1:p, col = mypalette(p), pch = 16, cex = 5)
# dev.off()


targets$colors <- mypalette(nlevels(targets$groups))[targets$groups]

write.table(targets, file = file.path("metadata", "targets.xls"), quote = FALSE, sep = "\t", row.names = FALSE)

########### read in targets

targets <- read.table(file.path("metadata", "targets.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)

targets


###########################################################################
#### import cel files
###########################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("pd.mogene.2.0.st")

library(oligo)
library(pd.mogene.2.0.st)

ff <- as.character(targets$FileName)

x <- oligo::read.celfiles(filenames = ff) ## GeneFeatureSet

pdf("boxplot.pdf")
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
boxplot(x, las = 2, col = targets$colors, names = targets$labels, las=2)
dev.off()

pdf("hist.pdf")
par(mar = c(5, 4, 4, 2) + 0.1)
hist(x, col = targets$colors, lwd = 2)
legend("topright", legend = targets$labels, col =  targets$colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()


###########################################################################
### PLM normalization; create images of chips, NUSE and RLE plots
###########################################################################

fitplm <- oligo::fitProbeLevelModel(x)


pdf("NUSE_fitplm.pdf")
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
oligo::NUSE(fitplm, col = targets$colors, names = targets$labels, las=2)
dev.off()

pdf("RLE_fitplm.pdf")
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
oligo::RLE(fitplm, col = targets$colors, names = targets$labels, las=2)
dev.off()




###########################################################################
### Normalization with RMA
###########################################################################


eset <- oligo::rma(x) ## Is the expression in log2 scale? ## ExpressionSet


pdf("boxplot_norm.pdf")
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
boxplot(eset, las = 2, col = targets$colors, names = targets$labels)
dev.off()

pdf("hist_norm.pdf")
par(mar = c(5, 4, 4, 2) + 0.1)
hist(eset, col = targets$colors, lwd = 2)
legend("topright", legend = targets$labels, col =  targets$colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()


labels <- paste0(ifelse(is.na(targets$ctrlRep), "", paste0(targets$ctrlRep, " ")), targets$labels)

pdf("MDS_norm.pdf")

### plot MDS for all samples
mds.dendo <- mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
### plot only points
plot(mds$x, mds$y, pch = 16, cex = 1.5, col = targets$colors)

### zoom in
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2, xlim = c(-1.5, 1.5), ylim = c(-1, 1))

### plot MDS for all samples except HeLa and whole bone marrow controls 
ex <- !targets$labels %in% c("control_wholeBoneMarrow", "control_HeLa")
mds <- plotMDS(eset[, ex], top=1000, col = targets$colors[ex], labels = labels[ex], cex = 1)

### plot only points
plot(mds$x, mds$y, pch = 16, cex = 1.5, col = targets$colors[ex] )

dev.off()




library(ggplot2)
library(ggdendro)

pdf("dendo_norm.pdf")

d <- mds.dendo$distance.matrix
rownames(d) <- labels
hc <- hclust(as.dist(d), method = "complete")
ggdendrogram(hc, theme_dendro = TRUE) + theme(axis.text.x = element_text(colour = targets$colors[hc$order], size = 13)) 

dev.off()



###########################################################################
####### NetAffx Annotation with getNetAffx()
###########################################################################


infoNetAffx <- pData(getNetAffx(eset, "transcript"))
head(infoNetAffx)



# apply(infoNetAffx, 2, function(cat){sum(is.na(cat))})
# 
# all(infoNetAffx$transcriptclusterid == infoNetAffx$probesetid)
# 
# sum(infoNetAffx$totalprobes)
# 
# ### check how many probesets have no annotation in fData and in infoNetAffx
# table(is.na(fData(eset)[,"ENTREZID"]))
# 
# table(is.na(fData(eset)[,"ENTREZID"]) & is.na(infoNetAffx$geneassignment))


###########################################################################
####### DO NOT RUN! Get annotation from NetAffx files from Affy website
### http://www.affymetrix.com/estore/browse/level_three_category_and_children.jsp?category=35868&categoryIdClicked=35868&expand=true&parent=35617
###########################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("AffyCompatible")

### Download the data 
library(AffyCompatible)

# password <- AffyCompatible:::acpassword

rsrc <- NetAffxResource(user = "gosia.nowicka@uzh.ch", password = "mockIP27", directory = "NetAffx")

availableArrays <- names(rsrc)
head(availableArrays)


availableArrays[grep("Mo", availableArrays)]


affxDescription(rsrc[["MoGene-2_0-st-v1"]])

annos <- rsrc[["MoGene-2_0-st-v1"]]
annos


anno <- affxAnnotation(annos)[[4]]
anno

fl <- readAnnotation(rsrc, annotation=anno, content=FALSE)




### Check what is in there
fl <- "NetAffx/MoGene-2_0-st-v1.na34.mm10.transcript.csv.zip"

conn <- unz(fl, "MoGene-2_0-st-v1.na34.mm10.transcript.csv")
# readLines(conn, n=20)

infoNetAffx2 <- read.table(conn, header=TRUE, sep=",", as.is = TRUE)
rownames(infoNetAffx2) <- infoNetAffx2$transcript_cluster_id

dim(infoNetAffx2)

apply(infoNetAffx2, 2, function(cat){sum(cat == "---")})



#### compare infoNetAffx2 with infoNetAffx
# all(infoNetAffx2$transcript_cluster_id == infoNetAffx2$probeset_id)
# colnames(infoNetAffx) <- colnames(infoNetAffx2)
# 
# probesetID <- "17457722" ## probe set with no ENTREZID
# infoNetAffx2[probesetID,]
# infoNetAffx[probesetID,]
# 
# infoNetAffx2[probesetID,] == infoNetAffx[probesetID,]
# 
# infoNetAffx2[probesetID, "mrna_assignment"]
# infoNetAffx[probesetID, "mrna_assignment"]


probesetID <- "17457722" ## probe set with no ENTREZID

infoNetAffx2[probesetID, "gene_assignment"]


geneAssi <- strsplit(infoNetAffx2$gene_assignment, " /// ")





###########################################################################
####### remove control probes == keep main probes
###########################################################################


# source("http://bioconductor.org/biocLite.R")
# biocLite("affycoretools")


# library(affycoretools)
# eset.main <- affycoretools::getMainProbes(eset) ### gives different results



table(infoNetAffx$category, useNA = "always")

all(featureNames(eset) == rownames(infoNetAffx))

keepMAIN <- infoNetAffx$category == "main"

eset.main <- eset[keepMAIN, ]


###########################################################################
####### Keep probes from chr1-chr19, Y, X
###########################################################################


table(infoNetAffx$seqname, useNA = "always")

keepCHR <- featureNames(eset.main) %in% rownames(infoNetAffx)[which(infoNetAffx$seqname %in% paste0("chr", c(1:19, "Y", "X")), useNames = TRUE)]

table(keepCHR)

eset.main <- eset.main[keepCHR, ]


###########################################################################
####### DO NOT USE THIS ONE. Annotation from mogene20sttranscriptcluster - has many entrez IDs missing
###########################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("mogene20sttranscriptcluster.db")


expr <- data.frame(exprs(eset))

library(mogene20sttranscriptcluster.db)

### Display all mappings
mogene20sttranscriptcluster()

# I way
annot <- data.frame(SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=" /// "), ENTREZID=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=" /// "), stringsAsFactors = FALSE)
colnames(annot) <- c("GeneSymbol_mogene20", "EntrezGeneID_mogene20")
annot[annot == "NA"] <- "---"

annot.mogene20 <- annot

annot.mogene20 <- annot.mogene20[featureNames(eset.main), ]

# # II way
# probes.ALL=row.names(expr)
# SYMBOL = unlist(mget(probes.ALL, mogene20sttranscriptclusterSYMBOL))
# ENTREZID = unlist(mget(probes.ALL, mogene20sttranscriptclusterENTREZID))
# 
# 
# ### check if it returns always one values - YES
# mg <- mget(probes.ALL, mogene20sttranscriptclusterENTREZID)
# table(sapply(mg, length))




# # IV way
# probes.ALL=row.names(expr)
# SYMBOLb = sapply(mget(probes.ALL, mogene20sttranscriptclusterSYMBOL), paste, collapse=", ")
# ENTREZIDb = sapply(mget(probes.ALL, mogene20sttranscriptclusterENTREZID), paste, collapse=", ")
# 
# all(SYMBOL == SYMBOLb)


# # III way
# library(annotate)
# probes.ALL <- featureNames(eset)
# SYMBOL <- getSYMBOL(probes.ALL,"mogene20sttranscriptcluster.db")

# annot <- data.frame(SYMBOL = SYMBOL, ENTREZID = ENTREZID , stringsAsFactors = FALSE)
# 
# fData(eset) <- annot 
# 
# 
# table(is.na(annot$ENTREZID))
# table(is.na(annot$SYMBOL))
# 
# 
# eset.org <- eset


###########################################################################
####### Get annotation from formated files from Affy website (NetAffx Analysis Center)
### http://www.affymetrix.com/analysis/index.affx
###########################################################################

library(limma)

###### files that are formated for easy load 


annof <- list.files("NetAffx", pattern = ".tsv", full.names = TRUE)
annof

### does not work
# anno_list <- read.table(annof[2], header = TRUE, sep = "\t", as.is = TRUE)



############## use public_database_references

# allLines <- readLines(annof[3], n=-1)
# 
# pdr <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
# colnames(pdr) <- gsub(pattern = " ", replacement = "" ,pdr[1,])
# pdr <- pdr[-1,]
# rownames(pdr) <- pdr$TranscriptClusterID
# 
# head(pdr)
# dim(pdr)
# 
# colnames(pdr)
# 
# table(pdr$EntrezGeneID == "---")
# table(pdr$GeneSymbol == "---")
# table(pdr$TranscriptID == "---")
# table(pdr$GOID == "---")
# 
# table(pdr[ featureNames(eset.main) ,"GOID"]== "---")
# table(pdr[ featureNames(eset.main) ,"TranscriptID"]== "---")
# table(pdr[ featureNames(eset.main) ,"GeneSymbol"] == "---")
# 
# probesetID <- "17299972" ## probe set with no ENTREZID
# 
# pdr[probesetID,]
# 
# pdr[probesetID, "GeneSymbol"]
# 
# infoNetAffx2[probesetID,]



############## use gene_list

allLines <- readLines(annof[2], n=-1)

gl <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(gl) <- gsub(pattern = " ", replacement = "" ,gl[1,])
gl <- gl[-1,]
rownames(gl) <- gl$TranscriptClusterID

colnames(gl)

### check for how many probe sets there is GO 

head(gl$GODescription)
table(gl$GODescription == "---")


# dim(gl)

table(gl$EntrezGeneID == "---")
table(gl$GeneSymbol == "---")
table(gl$GeneTitle == "---")


table(gl[ featureNames(eset.main) ,"GeneSymbol"] == "---")
table(gl[ featureNames(eset.main) ,"GeneTitle"] == "---")


### list of probe sets with no annotation in the end 
noAnnot <- featureNames(eset.main)[gl[ featureNames(eset.main) ,"GeneSymbol"] == "---"]

# probesetID <- "17422859" ## probe with ENTREZID: Tnfrsf4 22163
# infoNetAffx2[probesetID, 1:9]
# gl[probesetID,]




annot <- gl[ featureNames(eset.main) ,c("GeneSymbol", "EntrezGeneID", "GeneTitle")]
fData(eset.main) <- annot 



### compare annot with annot.mogene20 - weird thing for some probe sets the info is different... But what is in annot agrees with infoNetAffx2.

table(annot.mogene20$GeneSymbol_mogene20 == "---")
table(annot$GeneSymbol == "---")

table(annot.mogene20$EntrezGeneID_mogene20 == "---")
table(annot$EntrezGeneID == "---")

head(annot.mogene20)
head(annot)

infoNetAffx2["17210883", "gene_assignment"]
infoNetAffx["17210883", "geneassignment"]


infoNetAffx2["17210869", "gene_assignment"]
infoNetAffx["17210869", "geneassignment"]


infoNetAffx2["17210883", "mrna_assignment"]
infoNetAffx["17210883", "mrnaassignment"]



###########################################################################
####### Get ENSEMBL annotation using biomaRt
###########################################################################


library(biomaRt)


mart <- useMart("ensembl")
listDatasets(mart)


mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

attr <- listAttributes(mart)

attr[grep("affy", attr$name),]

# listFilters(mart)

genes <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "description","affy_mogene_2_1_st_v1"), filters="affy_mogene_2_1_st_v1", values=featureNames(eset.main), mart=mart)

dim(genes)
head(genes)

### clean the description
genes$description <- strsplit2(genes$description, " \\[Source")[, 1]


### some features have multiple ensembl annotation 
length(unique(genes$affy_mogene_2_1_st_v1))

probesetID <- "17457722" ## probe set with no ENTREZID
genes[genes$affy_mogene_2_1_st_v1 == probesetID, ]
gl[probesetID,]

probesetID <- "17422859" ## probe with ENTREZID: Tnfrsf4 22163
genes[genes$affy_mogene_2_1_st_v1 == probesetID, ]
gl[probesetID,]


### check what are the extra annotations that I get with Ensembl
noAnnotMart <- genes[genes$affy_mogene_2_1_st_v1 %in% noAnnot, ]
head(noAnnotMart)
## most of them are the predicted genes
table(grepl("predicted", noAnnotMart$description))
head(noAnnotMart[!grepl("predicted", noAnnotMart$description), ])
## for predicted genes the gene symbol starts with "Gm"
noAnnotMart[grepl("predicted", noAnnotMart$description), "external_gene_name" ]



### Merge the info about multiple genes into one string
library(plyr)

genes.m <- ddply(genes, "affy_mogene_2_1_st_v1", summarize, GeneSymbol_Ensembl = paste0(external_gene_name, collapse = " /// "), GeneTitle_Ensembl = paste0(description, collapse = " /// "), EnsemblGeneID = paste0(ensembl_gene_id, collapse = " /// "))
h(genes.m)

probesetID <- "17422859" ## probe with ENTREZID: Tnfrsf4 22163
genes.m[genes.m$affy_mogene_2_1_st_v1 == probesetID, ]



dim(annot)

dim(genes.m)


annot.m <- merge(annot, genes.m, by.x = 0, by.y = "affy_mogene_2_1_st_v1", all.x = TRUE, sort = FALSE, )
colnames(annot.m)[1] <- "ProbesetID"
rownames(annot.m) <- annot.m[,"ProbesetID"]

annot.m[is.na(annot.m)] <- "---"



### some checks
table(annot.m$GeneSymbol == annot.m$GeneSymbol_Ensembl)

table(annot.m$GeneSymbol == "---", !annot.m$GeneSymbol_Ensembl == "---")


head(annot.m[annot.m$GeneSymbol == "---" & !annot.m$GeneSymbol_Ensembl == "---", ])

extraAnnot <- !grepl("Gm",annot.m[, "GeneSymbol_Ensembl"]) & annot.m$GeneSymbol == "---" & !annot.m$GeneSymbol_Ensembl == "---"

table(extraAnnot)

annot.m[extraAnnot, c("GeneSymbol_Ensembl", "GeneTitle_Ensembl" )]


write.table(annot.m, "Annotation.xls", quote = FALSE, sep = "\t", row.names = FALSE)

fData(eset.main) <- annot.m



###########################################################################
####### Get probe info - probe 2 transcript cluster match
###########################################################################


### get probe 2 transcript match
probeInfo <- oligo::getProbeInfo(x, field = c('fid', 'fsetid', 'level', 'type', 'transcript_cluster_id'), probeType = "pm", target='core')

head(probeInfo)

table(probeInfo$type, useNA = "always")


setequal(featureNames(eset), unique(probeInfo$transcript_cluster_id))



###########################################################################
### Get GC content per probe
###########################################################################
# probe with higher GC content will have higher background

# pms <- oligo::pm(x, target='core')

pmSeq <- oligo::pmSequence(x, target='core')

library(Biostrings)
gcCont <- letterFrequency(pmSeq, letters='CG')[,1]

table(gcCont)

probeInfo$gcCont <- gcCont



###########################################################################
####### Filtering probes with low expression 
###########################################################################

#### using background information from antigenomic probesets

# library(genefilter)
# 
# tblNA <- table(infoNetAffx$category, useNA = "always")
# 
# antigm <- infoNetAffx[infoNetAffx$category == "control->bgp->antigenomic", "probesetid"]
# 
# bgExpr <- exprs(eset)[as.character(antigm), ]
# 
# 
# # bkg <- apply(bgExpr, 2, quantile, probs=0.5)
# # minval <- max(bkg)
# # minval
# 
# 
# bkg <- apply(bgExpr, 2, mean)
# # bkg <- rowMeans( bgExpr )
# 
# minval <- mean(bkg)
# minval
# 
# keep <- genefilter(eset.main, filterfun(kOverA(3, minval)))
# table(keep)
# 
# eset.main <- eset.main[keep,]




#################### based on GC content

### Get the background expression levels for different GC ammount

antigm <- infoNetAffx[infoNetAffx$category == "control->bgp->antigenomic", "probesetid"]

bgExpr <- exprs(eset)[as.character(antigm), targets$labels != "control_HeLa"]
bgExpr

bgProbeInfo <- subset(probeInfo, probeInfo$type == "control->bgp->antigenomic")

head(bgProbeInfo)

### see how many probes are for each GC content
table(bgProbeInfo$gcCont)


library(plyr)
library(ggplot2)
library(reshape2)

bgTransInfo <- ddply(bgProbeInfo, "transcript_cluster_id", summarize,  gcCont=mean(gcCont))


bgdf <- data.frame(bgTransInfo, bgExpr)

bgdf.m <- melt(bgdf, id.vars = c("transcript_cluster_id", "gcCont"), variable.name = "Samples", value.name = "Expression")
head(bgdf.m)

bgdf.m$gcCont <- factor(bgdf.m$gcCont)

ggp.bg <- ggplot(data = bgdf.m, aes(x = gcCont, y = Expression)) +
  geom_boxplot(colour = "lightcoral") +
  theme_bw()

pdf("gc_boxplot.pdf")
print(ggp.bg)
dev.off()


expr <- exprs(eset.main)[,targets$labels != "control_HeLa"]

### Get the GC content for all the probe sets
transInfo <- ddply(probeInfo, "transcript_cluster_id", summarize,  gcCont=mean(gcCont))
rownames(transInfo) <- transInfo$transcript_cluster_id
transInfo <- transInfo[rownames(expr), ]
transInfo$gcCont <- round(transInfo$gcCont)

table(transInfo$gcCont)

df <- data.frame(transInfo, expr)

df.m <- melt(df, id.vars = c("transcript_cluster_id", "gcCont"), variable.name = "Samples", value.name = "Expression")
head(df.m)


df.m$Type <- "Main"
bgdf.m$Type <- "BGP"

df.all <- rbind(df.m, bgdf.m)
df.all$gcCont <- factor(df.all$gcCont, levels = 3:25)


ggp <- ggplot(data = df.all, aes(x = gcCont, y = Expression, fill = Type)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position="top")

pdf("gc_boxplot_main_and_bgp.pdf")
print(ggp)
dev.off()



library(matrixStats)
# ls("package:matrixStats")

bgTransInfo$MedianExpr <- rowMedians(bgExpr)
bgTransInfo$Q075Expr <- rowQuantiles(bgExpr, probs = 0.75)
bgTransInfo


# pdf("gc.pdf")
# plot(bgTransInfo$gcCont, bgTransInfo$MedianExpr, type = "p", xlab = "GC content", ylab = "Median log2 expression", pch = 16, col = "lightcoral", cex = 2)
# dev.off()



transInfo$minExpr <- factor(transInfo$gcCont, levels = bgTransInfo$gcCont)
levels(transInfo$minExpr) <- bgTransInfo$MedianExpr
# levels(transInfo$minExpr) <- bgTransInfo$Q075Expr
transInfo$minExpr <- as.numeric(as.character(transInfo$minExpr))
head(transInfo)



all(rownames(expr) == transInfo$transcript_cluster_id)


keepEXPR <- sapply(1:nrow(expr), function(tr){ sum(expr[tr, ] > transInfo$minExpr[tr]) >= ncol(expr) } )

table(keepEXPR)


eset.main <- eset.main[keepEXPR, ]


eset.main


###########################################################################
#### Comparison 1: leukemia VS. controls
#### fitting model for all data
###########################################################################

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset.main <- eset.main[ ,ord]

# sampleNames(eset.main)
# targets$FileName

### keep only leukemia and control CD4+, CD4+CD8+ and CD8+ samples
samples2keep <- grepl("leukemia|control_CD", targets$labels)

treatments <- data.frame(Treatment = as.character(targets$groups[samples2keep]))
treatments

eset <- eset.main[, samples2keep] 

design <- model.matrix(~ 0 + Treatment, data=treatments)
rownames(design) <- targets$labels[samples2keep]
design


fit <- lmFit(eset, design)

contrasts <- cbind(CtrlCD4 = c(-1, 0, 0, 1), CtrlCD4CD8 = c(0, -1, 0, 1), CtrlCD8 = c(0, 0, -1, 1))
contrasts


fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2, trend = TRUE)


pdf("plotSA_trend.pdf")
plotSA(fit2)
dev.off()


# fit2 <- contrasts.fit(fit, contrasts)
# 
# fit2 <- eBayes(fit2)
# 
# pdf("plotSA.pdf")
# plotSA(fit2)
# dev.off()


results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=0)
summary(results)

## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)


write.fit(fit2, results = results, file = "results.xls", adjust = "BH")


pdf("vennDiagram.pdf")
vennDiagram(results,include=c("up", "down"))
dev.off()



pdf("plotMA.pdf")
limma::plotMA(fit2, coef = "CtrlCD4", status = results[, "CtrlCD4"])
limma::plotMA(fit2, coef = "CtrlCD4CD8", status = results[, "CtrlCD4CD8"])
limma::plotMA(fit2, coef = "CtrlCD8", status = results[, "CtrlCD8"])
dev.off()


### volcano plots
library(ggplot2)


coefs <- c("CtrlCD4", "CtrlCD4CD8", "CtrlCD8")

pdf("volcanoplot.pdf")

for(i in 1:length(coefs)){
coef <- coefs[i] 
table <- topTable(fit2, coef=coef, n=Inf)
table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
gg1 <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) + theme_bw() +ggtitle(coef) +
  theme(legend.position = "none") +
  xlab("log2 fold change") + ylab("-log10 p-value")

print(gg1)

}

dev.off()



### histograms of p-values and adjusted p-values
pdf(paste0("hist_pvs.pdf"))
for(i in 1:length(coefs)){
  
  coef <- coefs[i]  
  table <- topTable(fit2, coef=coef, n=Inf)
  hist(table$P.Value, breaks = 100, main = coef, xlab = "P-values")
  #hist(table$adj.P.Val, breaks = 100, main = coef, xlab = "Adjusted p-values")
  
}
dev.off()




### plot expression of top sign. genes/probesets
library(ggplot2)
library(reshape2)

topn <- 20
exp <- exprs(eset)
xs <- 1:ncol(exp)

for(i in 1:length(coefs)){
  # i = 1
  coef <- coefs[i]
  
  tt <- topTable(fit2, coef=coef, n=Inf, p.value=0.05, lfc=1)
  write.table(tt, paste0("topTable_",coef,".xls"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  print(coef)
  print(head(tt[, c("GeneSymbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")], topn))
  
  topp <- rownames(tt)[1:topn]

  pdf(paste0("topExpression_",coef,".pdf"))
  par(mfrow=c(2,2))
  
  for(i in 1:topn){
    
    plot(xs,exp[topp[i], ], xaxt = "n", ylab = "log2 Expression", xlab = "", pch = 16, cex = 2, col = targets$colors[samples2keep], main = paste0(topp[i]), las = 2)
    axis(side=1, at=xs, labels=NULL, las=2)
     
  }
  
  dev.off()
  

  df <- data.frame(Gene = topp, exp[topp,])
  df.m <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")
  ### keep order of genes as in tt
  df.m$Gene <- factor(df.m$Gene, levels = topp)
  ### add Entrez ID to the facet labels
  lab.fct <- paste0(topp, "\n", strsplit2(tt[topp, "GeneSymbol"], " /// ")[,1])
  levels(df.m$Gene) <- lab.fct
  
  ggp <- ggplot(df.m, aes(x = Sample, y = Expression)) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 10), plot.title = element_text(size = 16), strip.text.x = element_text(size = 10)) +
    scale_x_discrete(labels=targets$groups[samples2keep]) +
    labs(title = coef, y = "Log2 expression") +
    geom_bar(stat = "identity", colour = targets$colors[samples2keep], fill = targets$colors[samples2keep]) +
    facet_wrap(~ Gene, scales="free_y", ncol=4) 
  
  pdf(paste0("topExpressionBarPlot_",coef,".pdf"), 11, 11)
  print(ggp)    
  dev.off()
  
}



###########################################################################
#### GO analysis
###########################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")

library(topGO)
library(Rgraphviz)


affyLib <- "mogene20sttranscriptcluster.db"

coefs <- c("CtrlCD4", "CtrlCD4CD8", "CtrlCD8")
coef <- coefs[1]

tt <- topTable(fit2, coef=coef, n=Inf)
geneList <- rep(0, nrow(tt))
names(geneList) <- rownames(tt)
geneList[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1] <- 1
table(geneList)

### Function used to create new topGOdata object
fun.gene.sel <- function(geneList) {
  return(geneList <- ifelse(geneList==0, FALSE, TRUE))
}


allRes.Fisher.merged <- NULL
allRes.Fisher.elim.merged <- NULL


for(go in c("BP","MF","CC")){
  # go = "BP"
  
  sampleGOdata <- new("topGOdata", description = paste0("Simple session for ", coef), ontology = go, allGenes = geneList, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
  
  print(sampleGOdata)
  
  
  # Fisher's exact test which is based on gene counts, and a Kolmogorov-Smirnov like test which computes enrichment based on gene scores
  # For the method classic each GO category is tested independently
  
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultFisher.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
  
  
  pValues.Fisher <- score(resultFisher)
  topNodes.Fisher <- length(pValues.Fisher)
  
  pValues.Fisher.elim <- score(resultFisher.elim)
  topNodes.Fisher.elim <- length(pValues.Fisher.elim)
  
  
  allRes.Fisher <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = topNodes.Fisher) 
  
  allRes.Fisher$GO <- go   
  allRes.Fisher.merged <- rbind(allRes.Fisher.merged, allRes.Fisher)
  
  
  allRes.Fisher.elim <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = topNodes.Fisher.elim)      
  
  allRes.Fisher.elim$GO <- go 
  allRes.Fisher.elim.merged <- rbind(allRes.Fisher.elim.merged, allRes.Fisher.elim)
  
  
  
  pdf(paste("GO_", go, ".pdf", sep=""))
  showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
  showSigOfNodes(sampleGOdata, score(resultFisher.elim), firstSigNodes = 5, useInfo = 'all')
  dev.off()
  
  
}


allRes.Fisher.merged$classicFisher.adj <- p.adjust(allRes.Fisher.merged$classicFisher, method = "BH")
allRes.Fisher.merged <- allRes.Fisher.merged[order(allRes.Fisher.merged$classicFisher.adj, decreasing = FALSE), ]

head(allRes.Fisher.merged, 20)


write.table(allRes.Fisher.merged, paste("GO_Fisher_",coef,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)


allRes.Fisher.elim.merged$elimFisher.adj <- p.adjust(allRes.Fisher.elim.merged$elimFisher, method = "BH")
allRes.Fisher.elim.merged <- allRes.Fisher.elim.merged[order(allRes.Fisher.elim.merged$elimFisher.adj, decreasing = FALSE), ]

head(allRes.Fisher.elim.merged, 20)

write.table(allRes.Fisher.elim.merged, paste("GO_Fisher_elim_",coef,".xls", sep=""), sep="\t", row.names=F)





###########################################################################
#### Gene set enrichment analysis with C7  Immunologic genes sets
###########################################################################

# gene sets from MSigDB with ENTREZ IDs
load("MSigDB_v4_0/mouse_c7_v4.rdata")

Mm.c7[1]
mysets <- Mm.c7
length(mysets)

table(sapply(mysets, length))


### Create an Index for camera
annot <- fData(eset.main)
table(annot$EntrezGeneID == "---")

### Too slow
# EntrezGeneID <- strsplit(annot$EntrezGeneID, " /// ")
# Index <- lapply(mysets, function(ms){sapply(EntrezGeneID, function(eg){any(eg %in% ms)})})


EntrezGeneID <- strsplit2(annot$EntrezGeneID, " /// ")

nrow = nrow(EntrezGeneID)
ncol = ncol(EntrezGeneID)

Index <- lapply(mysets, function(ms){  
  eg <- matrix(EntrezGeneID %in% ms, nrow = nrow, ncol = ncol, byrow = FALSE)
  rowSums(eg) > 0 
  })


# ms <- mysets[[4]]
# eg <- matrix(EntrezGeneID %in% ms, nrow = nrow(EntrezGeneID), ncol = ncol(EntrezGeneID), byrow = FALSE)
# apply(eg, 2, sum)
# table(rowSums(eg) > 0)
# 
# table(Index[[4]])


# ms <- c(1, 2, 3)
# eg <- c(2, 4)
# any(eg %in% ms)


### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset.main <- eset.main[ ,ord]

### keep only leukemia and control CD4+, CD4+CD8+ and CD8+ samples
samples2keep <- grepl("leukemia|control_CD", targets$labels)

treatments <- data.frame(Treatment = as.character(targets$groups[samples2keep]))
treatments

eset <- eset.main[, samples2keep] 

design <- model.matrix(~ 0 + Treatment, data=treatments)
rownames(design) <- targets$labels[samples2keep]
design


contrasts <- cbind(CtrlCD4 = c(-1, 0, 0, 1), CtrlCD4CD8 = c(0, -1, 0, 1), CtrlCD8 = c(0, 0, -1, 1))
contrasts



gsea <- camera(y = eset, index=Index, design=design, contrast=contrasts[,"CtrlCD4"], trend.var=FALSE)

head(gsea)

table(gsea$FDR < 0.05)


gsea <- camera(y = eset, index=Index, design=design, contrast=contrasts[,"CtrlCD4CD8"], trend.var=TRUE)

head(gsea)

table(gsea$FDR < 0.05)



gsea <- camera(y = eset, index=Index, design=design, contrast=contrasts[,"CtrlCD8"], trend.var=TRUE)

head(gsea)

table(gsea$FDR < 0.05)













