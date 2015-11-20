###########################################################################
# Created 16 Oct 2015
# BioC 3.0

# DE analysis of Affymetrix Mouse Gene 2.0 ST arrays (pd.mogene.2.0.st)
# Additional replicates 

# Update 27 Oct 2015 
# Add pre versus after treatment analysis 

###########################################################################

setwd("/home/Shared/data/array/microarrays_edwin")

path_plots <- "Analysis_Mouse_Rebuttal/Plots/"
path_results <- "Analysis_Mouse_Rebuttal/"

dir.create(path_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(path_results, showWarnings = FALSE, recursive = TRUE)


###########################################################################
# create targets table with information about samples from Micro_array_sample_list.txt
###########################################################################

library(limma)


metadata <- read.table("metadata/Micro_array_sample_list.txt", sep = "\t")

targets <- data.frame(metadata , FileName = list.files("CEL/", pattern="IA201502" ,full.names = TRUE))

colnames(targets) <- c("Experiment", "SampleNr", "CellType", "FileName")

targets$ExperimentShort <- targets$Experiment
levels(targets$ExperimentShort) <- c("Bone marrow control"="control", "Kit control"="control", "Leukemia"="leukemia", "Leukemia after treatment" = "afterTreatment", "T cell control" = "control", "Thymocyte control" = "control")

targets$CellTypeShort <- targets$CellType
levels(targets$CellTypeShort) <- c("907" = "907", "907 - Post Dex" = "907", "B2M10" = "B2M10", "B2M2" = "B2M2", "B2M3" = "B2M3", "B2M3 Post dex" = "B2M3", "B3M3" = "B3M3", "B3M30" = "B3M30", "CD4 T cells spleen 1" = "CD4", "CD4 T cells spleen 2" = "CD4", "CD4 T cells spleen 3" = "CD4", "CD4+8+ DP Thymocytes 1" = "CD4+8+", "CD4+8+ DP Thymocytes 2" = "CD4+8+", "CD4+8+ DP Thymocytes 3" = "CD4+8+", "CD8 T cells spleen 1" = "CD8", "CD8 T cells spleen 2" = "CD8", "CD8 T cells spleen 3" = "CD8", "HeLa control" = "HeLa", "Whole bone marrow 1" = "wholeBoneMarrow", "Whole bone marrow 2" = "wholeBoneMarrow", "Whole bone marrow 3" = "wholeBoneMarrow")

targets$labels <- factor(paste(targets$ExperimentShort, targets$CellTypeShort, sep="_" ))

targets$groups <- targets$labels
levels(targets$groups)[grep(pattern = "leukemia", levels(targets$groups))] <- "leukemia"
levels(targets$groups)[grep(pattern = "afterTreatment", levels(targets$groups))] <- "afterTreatment"

targets$ctrlRep <- c(rep("", 12), rep(1:3, rep(4, 3)))

nlevels(targets$groups)
levels(targets$groups)


cbPalette <- c("#D55E00", "#F0E442","#56B4E9", "#009E73", "#0072B2","#CC79A7", "#999999")

pdf("Colors.pdf", width = 15)
barplot(rep(1, 7), col = cbPalette, names.arg = levels(targets$groups))
dev.off()

targets$colors <- cbPalette[targets$groups]

write.table(targets, file = file.path("metadata", "targets.xls"), quote = FALSE, sep = "\t", row.names = FALSE)

###########################################################################
# add samples from Micro_array_sample_list_rebuttal.xls
###########################################################################


targets_batch1 <- read.table(file.path("metadata", "targets.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)
targets_batch1$batch <- 1

targets_batch2 <- read.table(file.path("metadata", "Micro_array_sample_list_rebuttal.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)
targets_batch2$batch <- c(rep(2, 5), 1)

colors <- unique(targets_batch1[, c("groups", "colors")])
targets_batch2$colors <- colors$colors[match(targets_batch2$groups, colors$groups)] 



targets <- rbind(targets_batch1, targets_batch2)

write.table(targets, file = file.path("metadata", "targets_all.xls"), quote = FALSE, sep = "\t", row.names = FALSE)


###########################################################################
# read in all targets
###########################################################################


targets_org <- targets <- read.table(file.path("metadata", "targets_all.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)


###########################################################################
#### import cel files
###########################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("pd.mogene.2.0.st")

library(oligo)
library(pd.mogene.2.0.st)

ff <- as.character(targets$FileName)

x <- oligo::read.celfiles(filenames = ff) ## GeneFeatureSet

pdf(paste0(path_plots, "boxplot.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
boxplot(x, las = 2, col = targets$colors, names = targets$labels, las=2)
dev.off()

pdf(paste0(path_plots, "hist.pdf"))
par(mar = c(5, 4, 4, 2) + 0.1)
hist(x, col = targets$colors, lwd = 2)
legend("topright", legend = targets$labels, col =  targets$colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()


###########################################################################
### PLM normalization; create images of chips, NUSE and RLE plots
###########################################################################

fitplm <- oligo::fitProbeLevelModel(x)


pdf(paste0(path_plots,"NUSE_fitplm.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
oligo::NUSE(fitplm, col = targets$colors, names = targets$labels, las=2)
dev.off()


pdf(paste0(path_plots, "RLE_fitplm.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
oligo::RLE(fitplm, col = targets$colors, names = targets$labels, las=2)
dev.off()




###########################################################################
### Normalization with RMA
###########################################################################


eset_org <- eset <- oligo::rma(x) ## Is the expression in log2 scale? ## ExpressionSet


pdf(paste0(path_plots, "boxplot_norm.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
boxplot(eset, las = 2, col = targets$colors, names = targets$labels)
dev.off()

pdf(paste0(path_plots, "hist_norm.pdf"))
par(mar = c(5, 4, 4, 2) + 0.1)
hist(eset, col = targets$colors, lwd = 2)
legend("topright", legend = targets$labels, col =  targets$colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()



###########################################################################
### MDS plots
###########################################################################

library(limma)


########## All samples

eset <- eset_org
targets <- targets_org

labels <- targets$groups


pdf(paste0(path_plots, "MDS_all.pdf"), width = 5, height = 5)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
dev.off()


legend <- unique(targets[, c("groups", "colors")])
min <- min(mds$x, mds$y)
max <- max(mds$x, mds$y)

pdf(paste0(path_plots, "MDS_all_points.pdf"), width = 5, height = 5)

plot(mds$x, mds$y, pch = targets$batch, col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(min, max), ylim = c(min, max))

# text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.3, cex = 0.3)

legend("topleft", legend = c(legend$groups, "batch1", "batch2"), pch = c(rep(16, nrow(legend)), 1, 2), col = c(legend$colors, 1, 1), cex = 0.8, bty = "n")

dev.off()



########## All samples with no hela 

keep_samps <- targets_org$CellTypeShort != "HeLa" 

eset <- eset_org[, keep_samps]
targets <- targets_org[keep_samps, ]

labels <- targets$groups


pdf(paste0(path_plots, "MDS_all_noHela.pdf"), width = 5, height = 5)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
dev.off()


legend <- unique(targets[, c("groups", "colors")])
min <- min(mds$x, mds$y)
max <- max(mds$x, mds$y)


pdf(paste0(path_plots, "MDS_all_noHela_points.pdf"), width = 5, height = 5)

plot(mds$x, mds$y, pch = targets$batch, col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(min, max), ylim = c(min, max))

# text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.3, cex = 0.3)

legend("bottomleft", legend = c(legend$groups, "batch1", "batch2"), pch = c(rep(16, nrow(legend)), 1, 2), col = c(legend$colors, 1, 1), cex = 0.8, bty = "n")

dev.off()

### zoom on leukemia and treatment samples

keep_samps <- targets$ExperimentShort != "control"

eset <- eset[, keep_samps]
targets <- targets[keep_samps, ]

legend <- unique(targets[, c("groups", "colors")])
min <- min(mds$x[keep_samps], mds$y[keep_samps])
max <- max(mds$x[keep_samps], mds$y[keep_samps])


pdf(paste0(path_plots, "MDS_all_noHela_points_zoom.pdf"), width = 5, height = 5)

plot(mds$x[keep_samps], mds$y[keep_samps], pch = targets$batch, col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(min, max), ylim = c(min, max))

text(mds$x[keep_samps], mds$y[keep_samps], labels = targets$CellTypeShort, pos = 3, offset = 0.4, cex = 0.5, col = targets$colors)

legend("bottomleft", legend = c(legend$groups, "batch1", "batch2"), pch = c(rep(16, nrow(legend)), 1, 2), col = c(legend$colors, 1, 1), cex = 0.8, bty = "n")

dev.off()




########## Only controls 

keep_samps <- targets_org$CellTypeShort != "HeLa" & targets_org$ExperimentShort == "control"

eset <- eset_org[, keep_samps]
targets <- targets_org[keep_samps, ]

labels <- targets$groups


pdf(paste0(path_plots, "MDS_controls.pdf"), width = 5, height = 5)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
dev.off()


legend <- unique(targets[, c("groups", "colors")])
min <- min(mds$x, mds$y)
max <- max(mds$x, mds$y)


pdf(paste0(path_plots, "MDS_controls_points.pdf"), width = 5, height = 5)

plot(mds$x, mds$y, pch = targets$batch, col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(min, max), ylim = c(min, max))

# text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.3, cex = 0.3)

legend("topright", legend = c(legend$groups), pch = c(rep(16, nrow(legend))), col = c(legend$colors), cex = 0.8, bty = "n")

dev.off()



########## Pre and after treatment

keep_samps <- targets_org$CellTypeShort != "HeLa" & targets_org$ExperimentShort != "control"

eset <- eset_org[, keep_samps]
targets <- targets_org[keep_samps, ]

labels <- targets$groups


pdf(paste0(path_plots, "MDS_treatment.pdf"), width = 5, height = 5)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
dev.off()


legend <- unique(targets[, c("groups", "colors")])
min <- min(mds$x, mds$y)
max <- max(mds$x, mds$y)


pdf(paste0(path_plots, "MDS_treatment_points.pdf"), width = 5, height = 5)

plot(mds$x, mds$y, pch = targets$batch, col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(min, max), ylim = c(min, max))

text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.4, cex = 0.5, col = targets$colors)

legend("bottomleft", legend = c(legend$groups, "batch1", "batch2"), pch = c(rep(16, nrow(legend)), 1, 2), col = c(legend$colors, 1, 1), cex = 0.8, bty = "n")


dev.off()





###########################################################################
####### Do NOT keep the HeLa sample for the rest of the analysis
###########################################################################



keepSAMPS <- targets_org$labels != "control_HeLa"

eset_org <- eset <- eset_org[, keepSAMPS]
targets_org <- targets <- targets_org[keepSAMPS, ]

save(targets_org, file = paste0(path_results, "targets_org.Rdata"))

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


# probesetID <- "17457722" ## probe set with no ENTREZID
# infoNetAffx2[probesetID, "gene_assignment"]
# geneAssi <- strsplit(infoNetAffx2$gene_assignment, " /// ")


###########################################################################
####### remove control probes == keep main probes
###########################################################################


# source("http://bioconductor.org/biocLite.R")
# biocLite("affycoretools")


# library(affycoretools)
# eset_main <- affycoretools::getMainProbes(eset) ### gives different results



table(infoNetAffx$category, useNA = "always")

all(featureNames(eset) == rownames(infoNetAffx))

keepMAIN <- infoNetAffx$category == "main"

eset_main <- eset[keepMAIN, ]


###########################################################################
####### Keep probes from chr1-chr19, Y, X
###########################################################################


table(infoNetAffx$seqname, useNA = "always")

keepCHR <- featureNames(eset_main) %in% rownames(infoNetAffx)[which(infoNetAffx$seqname %in% paste0("chr", c(1:19, "Y", "X")), useNames = TRUE)]

table(keepCHR)

eset_main <- eset_main[keepCHR, ]



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

annot_mergeogene20 <- annot

annot_mergeogene20 <- annot_mergeogene20[featureNames(eset_main), ]

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
# eset_org <- eset


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
# table(pdr[ featureNames(eset_main) ,"GOID"]== "---")
# table(pdr[ featureNames(eset_main) ,"TranscriptID"]== "---")
# table(pdr[ featureNames(eset_main) ,"GeneSymbol"] == "---")
# 
# probesetID <- "17299972" ## probe set with no ENTREZID
# 
# pdr[probesetID,]
# 
# pdr[probesetID, "GeneSymbol"]
# 
# infoNetAffx2[probesetID,]



############## use gene_list

allLines <- readLines(annof[grepl("gene_list", annof)], n=-1)

gl <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(gl) <- gsub(pattern = " ", replacement = "" ,gl[1,])
gl <- gl[-1,]
rownames(gl) <- gl$TranscriptClusterID

colnames(gl)



# ### check for how many probe sets there is GO 
# head(gl$GODescription)
# table(gl$GODescription == "---")
# 
# 
# # dim(gl)
# 
# table(gl$EntrezGeneID == "---")
# table(gl$GeneSymbol == "---")
# table(gl$GeneTitle == "---")
# 
# 
# table(gl[ featureNames(eset_main) ,"GeneSymbol"] == "---")
# table(gl[ featureNames(eset_main) ,"GeneTitle"] == "---")
# 
# 
# ### list of probe sets with no annotation in the end 
# noAnnot <- featureNames(eset_main)[gl[ featureNames(eset_main) ,"GeneSymbol"] == "---"]
# 
# # probesetID <- "17422859" ## probe with ENTREZID: Tnfrsf4 22163
# # infoNetAffx2[probesetID, 1:9]
# # gl[probesetID,]



annot <- gl[ featureNames(eset_main) ,c("GeneSymbol", "EntrezGeneID", "GeneTitle")]



# ### compare annot with annot_mergeogene20 - weird thing for some probe sets the info is different... But what is in annot agrees with infoNetAffx2.
# 
# table(annot_mergeogene20$GeneSymbol_mogene20 == "---")
# table(annot$GeneSymbol == "---")
# 
# table(annot_mergeogene20$EntrezGeneID_mogene20 == "---")
# table(annot$EntrezGeneID == "---")
# 
# head(annot_mergeogene20)
# head(annot)
# 
# infoNetAffx2["17210883", "gene_assignment"]
# infoNetAffx["17210883", "geneassignment"]
# 
# 
# infoNetAffx2["17210869", "gene_assignment"]
# infoNetAffx["17210869", "geneassignment"]
# 
# 
# infoNetAffx2["17210883", "mrna_assignment"]
# infoNetAffx["17210883", "mrnaassignment"]
# 
# 
# probeID <- "17532811" # Foxp3
# 
# annot_mergeogene20[probeID, ]
# annot[probeID, ]
# annot_merge[probeID, ]

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

genes <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "description","affy_mogene_2_1_st_v1"), filters="affy_mogene_2_1_st_v1", values=featureNames(eset_main), mart=mart)

dim(genes)
head(genes)

### clean the description
genes$description <- strsplit2(genes$description, " \\[Source")[, 1]



# ### Do some checks
# ### some features have multiple ensembl annotation 
# length(unique(genes$affy_mogene_2_1_st_v1))
# 
# probesetID <- "17457722" ## probe set with no ENTREZID
# genes[genes$affy_mogene_2_1_st_v1 == probesetID, ]
# gl[probesetID,]
# 
# probesetID <- "17422859" ## probe with ENTREZID: Tnfrsf4 22163
# genes[genes$affy_mogene_2_1_st_v1 == probesetID, ]
# gl[probesetID,]
# 
# 
# ### check what are the extra annotations that I get with Ensembl
# noAnnotMart <- genes[genes$affy_mogene_2_1_st_v1 %in% noAnnot, ]
# head(noAnnotMart)
# ## most of them are the predicted genes
# table(grepl("predicted", noAnnotMart$description))
# head(noAnnotMart[!grepl("predicted", noAnnotMart$description), ])
# ## for predicted genes the gene symbol starts with "Gm"
# noAnnotMart[grepl("predicted", noAnnotMart$description), "external_gene_name" ]



### Merge the info about multiple genes into one string
library(plyr)

genes_merge <- plyr::ddply(genes, "affy_mogene_2_1_st_v1", summarize, GeneSymbol_Ensembl = paste0(external_gene_name, collapse = " /// "), GeneTitle_Ensembl = paste0(description, collapse = " /// "), EnsemblGeneID = paste0(ensembl_gene_id, collapse = " /// "))
h(genes_merge)


# ### Do some checks
# probesetID <- "17422859" ## probe with ENTREZID: Tnfrsf4 22163
# genes_merge[genes_merge$affy_mogene_2_1_st_v1 == probesetID, ]
# dim(annot)
# 
# dim(genes_merge)


annot_merge <- merge(annot, genes_merge, by.x = 0, by.y = "affy_mogene_2_1_st_v1", all.x = TRUE, sort = FALSE)
colnames(annot_merge)[1] <- "ProbesetID"
rownames(annot_merge) <- annot_merge[,"ProbesetID"]

annot_merge[is.na(annot_merge)] <- "---"



# ### some checks
# table(annot_merge$GeneSymbol == annot_merge$GeneSymbol_Ensembl)
# 
# table(annot_merge$GeneSymbol == "---", !annot_merge$GeneSymbol_Ensembl == "---")
# 
# head(annot_merge[annot_merge$GeneSymbol == "---" & !annot_merge$GeneSymbol_Ensembl == "---", ])
# 
# extraAnnot <- !grepl("Gm",annot_merge[, "GeneSymbol_Ensembl"]) & annot_merge$GeneSymbol == "---" & !annot_merge$GeneSymbol_Ensembl == "---"
# 
# table(extraAnnot)
# 
# annot_merge[extraAnnot, c("GeneSymbol_Ensembl", "GeneTitle_Ensembl" )]



all(annot_merge$ProbesetID == featureNames(eset_main))


fData(eset_main) <- annot_merge[featureNames(eset_main), ]




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
# keep <- genefilter(eset_main, filterfun(kOverA(3, minval)))
# table(keep)
# 
# eset_main <- eset_main[keep,]




#################### based on GC content

### Get the background expression levels for different GC ammount

antigm <- infoNetAffx[infoNetAffx$category == "control->bgp->antigenomic", "probesetid"]

bgExpr <- exprs(eset)[as.character(antigm), ]
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

pdf(paste0(path_plots, "gc_boxplot.pdf"))
print(ggp.bg)
dev.off()


expr <- exprs(eset_main)

### Get the GC content for all the probe sets
transInfo <- ddply(probeInfo, "transcript_cluster_id", summarize,  gcCont=mean(gcCont))
rownames(transInfo) <- transInfo$transcript_cluster_id
transInfo <- transInfo[rownames(expr), ]
transInfo$gcCont <- round(transInfo$gcCont)
### see what is the average GC content for main probe sets
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

pdf(paste0(path_plots, "gc_boxplot_main_and_bgp.pdf"))
print(ggp)
dev.off()


################### set the threshold for each probe set
library(matrixStats)
# ls("package:matrixStats")

bgTransInfo$Q095Expr <- rowQuantiles(bgExpr, probs = 0.95)
bgTransInfo


# pdf(paste0(path_plots, "gc.pdf"))
# plot(bgTransInfo$gcCont, bgTransInfo$MedianExpr, type = "p", xlab = "GC content", ylab = "Median log2 expression", pch = 16, col = "lightcoral", cex = 2)
# dev.off()

transInfo$minExpr <- factor(transInfo$gcCont, levels = bgTransInfo$gcCont)
levels(transInfo$minExpr) <- bgTransInfo$Q095Expr
transInfo$minExpr <- as.numeric(as.character(transInfo$minExpr))
head(transInfo)

save(transInfo, file = paste0(path_results, "transInfo.Rdata"))


#### Filtering itself

all(rownames(expr) == transInfo$transcript_cluster_id)


keepEXPR <- sapply(1:nrow(expr), function(tr){ sum(expr[tr, ] > transInfo$minExpr[tr]) >= 3 } )

table(keepEXPR)


eset_main <- eset_main[keepEXPR, ]


eset_main_org <- eset_main

save(eset_main_org, file = paste0(path_results, "eset_main_org.Rdata"))


###########################################################################
##### Multiple plot function
###########################################################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


###########################################################################
#### Comparison 1: leukemia VS. controls
#### fitting model for all data
###########################################################################

library(oligo)
library(pd.mogene.2.0.st)
library(limma)

load(paste0(path_results, "eset_main_org.Rdata"))
load(paste0(path_results, "targets_org.Rdata"))


targets <- targets_org
eset_main <- eset_main_org

### keep only leukemia and control CD4+, CD4+CD8+ and CD8+ and bone marrow samples
samples2keep <- grepl("leukemia|control", targets$labels)

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset_main <- eset_main[ ,ord]

# all(sampleNames(eset_main) == strsplit2(targets$FileName, "//")[,2])


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups))
treatments

design <- model.matrix(~ 0 + Treatment, data=treatments)
rownames(design) <- targets$labels
design

fit <- lmFit(eset_main, design)

contrasts <- cbind(CtrlCD4 = c(-1, 0, 0, 0, 1), CtrlCD4CD8 = c(0, -1, 0, 0, 1), CtrlCD8 = c(0, 0, -1, 0, 1), CtrlBM = c(0, 0, 0, -1, 1)) # treatment - control
contrasts


fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2, trend = TRUE)


pdf(paste0(path_plots, "Comp1_plotSA_trend.pdf"))
plotSA(fit2)
dev.off()


## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)



colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp1_vennDiagram.pdf"))
vennDiagram(results,include=c("up", "down"), circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="both", circle.col=colours, counts.col=colours)
# vennDiagram(results,include="up", circle.col=colours, counts.col=colours)
# vennDiagram(results,include="down", circle.col=colours, counts.col=colours)
dev.off()



### save all results with nice order
coefs <- c("CtrlCD4", "CtrlCD4CD8", "CtrlCD8", "CtrlBM")

# resExpr <- round(exprs(eset_main), 2)
# colnames(resExpr) <- paste0(treatments$Treatment, "_", colnames(resExpr))

resExpr <- round(exprs(eset_main_org), 2)
colnames(resExpr) <- paste0(targets_org$labels, "_", colnames(resExpr))
resExpr <- resExpr[, order(colnames(resExpr))]

resCoeff <- fit2$coefficients
colnames(resCoeff) <- paste0(colnames(resCoeff), "_coeffs")
resT <- fit2$t
colnames(resT) <- paste0(colnames(resT), "_t")
resPValue <- fit2$p.value
colnames(resPValue) <- paste0(colnames(resPValue), "_PValues")
resPValueAdj <- apply(fit2$p.value, 2, p.adjust, method = "BH")
colnames(resPValueAdj) <- paste0(colnames(resPValueAdj), "_AdjPValues")
resGenes <- fit2$genes
resRes <- results
colnames(resRes) <- paste0(colnames(resRes), "_Results")

stats <- c("coeffs", "t", "PValues", "AdjPValues", "Results")
colOrder <- paste(rep(coefs, each = length(stats)), rep(stats, length(coefs)), sep = "_")


resDE <- data.frame(resCoeff, resT, resPValue, resPValueAdj, resRes)[, colOrder]

resAll <- cbind(resGenes, resDE,  resExpr)

write.table(resAll, file = paste0(path_results, "Comp1_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)





### plot MA
pdf(paste0(path_plots, "Comp1_plotMA.pdf"))

for(i in 1:length(coefs)){
  coef <- coefs[i]
  limma::plotMA(fit2, coef = coef, status = results[, coef], values = c(-1, 0, 1), col = c("red", "black", "green"), cex = c(0.7, 0.3, 0.7), main = coef)
  abline(0,0,col="blue")
}

dev.off()




### volcano plots
library(ggplot2)

coefs <- c("CtrlCD4", "CtrlCD4CD8", "CtrlCD8", "CtrlBM")

gg1 <- list()
for(i in 1:length(coefs)){
  coef <- coefs[i] 
  table <- topTable(fit2, coef=coef, n=Inf)
  table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
  gg1[[i]] <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) + theme_bw() +ggtitle(coef) +
    theme(legend.position = "none") +
    xlab("log2 fold change") + ylab("-log10 p-value")
  
}

pdf(paste0(path_plots, "Comp1_volcanoplot.pdf"))
print(multiplot(plotlist = gg1, cols=2))
dev.off()



### histograms of p-values and adjusted p-values
colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp1_hist_pvs.pdf"))
for(i in 1:length(coefs)){
  
  coef <- coefs[i]  
  table <- topTable(fit2, coef=coef, n=Inf)
  hist(table$P.Value, breaks = 100, main = coef, xlab = "P-values", col = colours[i])
  #hist(table$adj.P.Val, breaks = 100, main = coef, xlab = "Adjusted p-values")
  
}
dev.off()




# ### plot expression of top sign. genes/probesets
# library(ggplot2)
# library(reshape2)
# 
# topn <- 20
# expr <- exprs(eset_main)
# xs <- 1:ncol(expr)
# 
# for(i in 1:length(coefs)){
#   # i = 1
#   coef <- coefs[i]
#   print(coef)
#   
#   tt <- topTable(fit2, coef=coef, n=Inf, p.value=0.05, lfc=1)
#   # write.table(tt, paste0("Comp1_topTable_",coef,".xls"), quote = FALSE, sep = "\t", row.names = FALSE)
#   
#   ### in the report display only first gene symbol
#   GeneSymbol <- strsplit2(head(tt[,"GeneSymbol"], topn), " /// ")[,1]
#   GeneTitle <- paste0(substr(strsplit2(head(tt[,"GeneTitle"], topn), " /// ")[,1], 1, 30))
#   
#   print(data.frame(GeneSymbol = GeneSymbol, GeneTitle = GeneTitle , head(tt[, c("logFC", "AveExpr", "P.Value", "adj.P.Val")], topn)))
#   
#   topp <- rownames(tt)[1:topn]
#   
#   df <- data.frame(Gene = topp, expr[topp,])
#   df.m <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")
#   ### keep order of genes as in tt
#   df.m$Gene <- factor(df.m$Gene, levels = topp)
#   ### add Entrez ID to the facet labels
#   lab.fct <- paste0(topp, "\n", strsplit2(tt[topp, "GeneSymbol"], " /// ")[,1])
#   levels(df.m$Gene) <- lab.fct
#   
#   ggp <- ggplot(df.m, aes(x = Sample, y = Expression)) +  
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 10), plot.title = element_text(size = 16), strip.text.x = element_text(size = 10)) +
#     scale_x_discrete(labels=targets$groups) +
#     labs(title = coef, y = "Log2 expression") +
#     geom_bar(stat = "identity", colour = targets$colors, fill = targets$colors) +
#     facet_wrap(~ Gene, scales="free_y", ncol=4) 
#   
#   pdf(paste0(path_plots, "Comp1_topExpressionBarPlot_",coef,".pdf"), 11, 11)
#   print(ggp)    
#   dev.off()
#   
# }



###########################################################################
#### Gene set enrichment analysis with C5 - GO genes sets
###########################################################################

# gene sets from MSigDB with ENTREZ IDs
load("MSigDB_v4_0/mouse_c5_v4.rdata")

mysets <- Mm.c5
length(mysets)

### keep the sets of interest
intrset <- read.table("Gene_Sets/Interesting_gene_sets_C5.txt", header = FALSE, sep = ",")[, 1]
intrset

intrset <- gsub("-", " ", intrset)
intrset <- gsub(" ", "_", intrset)

intrset <- toupper(intrset)
length(intrset)

sum(names(mysets) %in% intrset)

mysets <- mysets[intrset]


# table(sapply(mysets, length))


### Create an Index for camera
annot <- fData(eset_main)
# table(annot$EntrezGeneID == "---")

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


IndexMx <- do.call(cbind, Index)
class(IndexMx) <- "numeric"
colnames(IndexMx) <- names(mysets)
IndexMx <- data.frame(ProbesetID = annot$ProbesetID, IndexMx)

resAll <- merge(resAll, IndexMx, by = "ProbesetID", sort = FALSE)

write.table(resAll, file = "Comp1_DE_results_AllPlus.xls", quote = FALSE, sep = "\t", row.names = FALSE)


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups))

design <- model.matrix(~ 0 + Treatment, data=treatments)
rownames(design) <- targets$labels
design

contrasts <- cbind(CtrlCD4 = c(-1, 0, 0, 0, 1), CtrlCD4CD8 = c(0, -1, 0, 0, 1), CtrlCD8 = c(0, 0, -1, 0, 1), CtrlBM = c(0, 0, 0, -1, 1)) # treatment - control
contrasts




### run CAMERA

gsea <- list()

coef <- "CtrlCD4"
gsea.tmp <- gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=FALSE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("NGenes","Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), NGenes = gsea[[coef]][,1], gsea[[coef]][,-1])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)



### using information from eBayes fitting: fit2

pdf(paste0("PLOTS/GS_barcodeplot_",coef,".pdf"))

topgs <- 1
gsn <-rownames(gsea[[coef]])[1:topgs] 
gss <- gsea.tmp[gsn, , drop = FALSE]

for(i in 1:length(topgs)){
  
  barcodeplot(statistics = as.numeric((fit2$t[,coef])), index = Index[[gsn[i]]], index2 = NULL, gene.weights = as.numeric((fit2$coefficients[, coef]))[Index[[gsn[i]]]], weights.label = "logFC", labels = c("Up","Down"), quantiles = c(-1,1), col.bars = NULL, worm = TRUE, span.worm=0.45, main = paste0(gsn[i], "\n", gss[i, "Direction"], ", FDR = ", sprintf("%.02e",gss[i, "FDR"])))
  
  barcodeplot(statistics = as.numeric((fit2$p.value[, coef])), index = Index[[gsn[i]]], index2 = NULL, gene.weights = as.numeric((fit2$coefficients[, coef]))[Index[[gsn[i]]]], weights.label = "logFC", labels = c("Not significant","Significant"), quantiles = c(0.05,1), col.bars = NULL, worm = TRUE, span.worm=0.45, main = paste0(gsn[i], "\n", gss[i, "Direction"], ", FDR = ", sprintf("%.02e",gss[i, "FDR"])))
  
  
}

dev.off()




coef <- "CtrlCD4CD8"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)


coef <- "CtrlCD8"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)



coef <- "CtrlBM"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)





### merge all results into one table
gseaAll <- merge(gsea[["CtrlCD4"]], gsea[["CtrlCD4CD8"]], by = "GeneSet", all = TRUE)
gseaAll <- merge(gseaAll, gsea[["CtrlCD8"]], by = "GeneSet", all = TRUE)
gseaAll <- merge(gseaAll, gsea[["CtrlBM"]], by = "GeneSet", all = TRUE)
write.table(gseaAll, paste("Comp1_GSEA_C5_All.xls", sep=""), sep="\t", row.names=F, quote = FALSE)





# http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm



###########################################################################
#### Gene set enrichment analysis with C7  Immunologic genes sets
###########################################################################

# gene sets from MSigDB with ENTREZ IDs
load("MSigDB_v4_0/mouse_c7_v4.rdata")

# Mm.c7[1]
mysets <- Mm.c7
length(mysets)

# table(sapply(mysets, length))


### Create an Index for camera
annot <- fData(eset_main)
# table(annot$EntrezGeneID == "---")

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


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups))

design <- model.matrix(~ 0 + Treatment, data=treatments)
rownames(design) <- targets$labels
design

contrasts <- cbind(CtrlCD4 = c(-1, 0, 0, 0, 1), CtrlCD4CD8 = c(0, -1, 0, 0, 1), CtrlCD8 = c(0, 0, -1, 0, 1), CtrlBM = c(0, 0, 0, -1, 1)) # treatment - control
contrasts




### run CAMERA

gsea <- list()

coef <- "CtrlCD4"
gsea.tmp <- gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=FALSE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("NGenes","Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), NGenes = gsea[[coef]][,1], gsea[[coef]][,-1])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)



### using information from eBayes fitting: fit2

pdf(paste0("PLOTS/GS_barcodeplot_",coef,".pdf"))

topgs <- 1
gsn <-rownames(gsea[[coef]])[1:topgs] 
gss <- gsea.tmp[gsn, , drop = FALSE]

for(i in 1:length(topgs)){
  
  barcodeplot(statistics = as.numeric((fit2$t[,coef])), index = Index[[gsn[i]]], index2 = NULL, gene.weights = as.numeric((fit2$coefficients[, coef]))[Index[[gsn[i]]]], weights.label = "logFC", labels = c("Up","Down"), quantiles = c(-1,1), col.bars = NULL, worm = TRUE, span.worm=0.45, main = paste0(gsn[i], "\n", gss[i, "Direction"], ", FDR = ", sprintf("%.02e",gss[i, "FDR"])))
  
  barcodeplot(statistics = as.numeric((fit2$p.value[, coef])), index = Index[[gsn[i]]], index2 = NULL, gene.weights = as.numeric((fit2$coefficients[, coef]))[Index[[gsn[i]]]], weights.label = "logFC", labels = c("Not significant","Significant"), quantiles = c(0.05,1), col.bars = NULL, worm = TRUE, span.worm=0.45, main = paste0(gsn[i], "\n", gss[i, "Direction"], ", FDR = ", sprintf("%.02e",gss[i, "FDR"])))
  
  
}

dev.off()




coef <- "CtrlCD4CD8"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)


coef <- "CtrlCD8"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)



coef <- "CtrlBM"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)





### merge all results into one table
gseaAll <- merge(gsea[["CtrlCD4"]], gsea[["CtrlCD4CD8"]], by = "GeneSet", all = TRUE)
gseaAll <- merge(gseaAll, gsea[["CtrlCD8"]], by = "GeneSet", all = TRUE)
gseaAll <- merge(gseaAll, gsea[["CtrlBM"]], by = "GeneSet", all = TRUE)


write.table(gseaAll, paste("Comp1_GSEA_C7_All.xls", sep=""), sep="\t", row.names=F, quote = FALSE)





# http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm




###########################################################################
#### Gene set enrichment analysis with Hallmark genes sets
###########################################################################

############### Create mouse_hallmark_v5.rdata object like on WEHI web

allLines <- readLines("MSigDB_v4_0/h.all.v5.0.entrez.gmt", n = -1)

humanSets <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)

namesHS <- humanSets[, 1]

Hu.hallmark <- apply(humanSets, 1, function(r){ 
  r <- r[-c(1,2)]
  r <- r[r != ""]
  return(as.numeric(r))
} )


names(Hu.hallmark) <- namesHS

### get the mouse human homology
hom <- read.table("MSigDB_v4_0/HOM_MouseHumanSequence.txt", header = TRUE, sep = "\t")

homM <- hom[hom$Common.Organism.Name == "mouse, laboratory", c("HomoloGene.ID", "EntrezGene.ID")]

homH <- hom[hom$Common.Organism.Name == "human", c("HomoloGene.ID", "EntrezGene.ID")]

homMatch <- merge(homH, homM, by = "HomoloGene.ID", sort = FALSE, all = TRUE) 

homMatch <- homMatch[!is.na(homMatch[, 2]) & !is.na(homMatch[, 3]), ]

# merge(data.frame(a = c(1, 1, 2), b = c(21, 23, 24)), data.frame(a = c(1, 1, 2, 2), b = c(31, 32, 33, 34)) , by=1, sort = FALSE, all = TRUE)


Mm.hallmark <- lapply(Hu.hallmark, function(gs){
  
  unique(homMatch[homMatch[, 2] %in% gs, 3])
  
})


save(Mm.hallmark, file = "MSigDB_v4_0/mouse_hallmark_v5.rdata")


############### Create mouse_hallmark_v5.rdata object like on WEHI web


# gene sets from MSigDB with ENTREZ IDs
load("MSigDB_v4_0/mouse_hallmark_v5.rdata")

mysets <- Mm.hallmark
length(mysets)

### keep the sets of interest
intrset <- read.table("Gene_Sets/Interesting_gene_sets_Hallmark.txt", header = FALSE, sep = ",", as.is = TRUE)[, 1]
intrset

intrset <- gsub("-", " ", intrset)
intrset <- gsub(" ", "_", intrset)

intrset <- paste0("HALLMARK_",toupper(intrset))
length(intrset)

sum(names(mysets) %in% intrset)

# intrset[!intrset %in% names(mysets)]


mysets <- mysets[intrset]


# table(sapply(mysets, length))


### Create an Index for camera
annot <- fData(eset_main)
# table(annot$EntrezGeneID == "---")

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


IndexMx <- do.call(cbind, Index)
class(IndexMx) <- "numeric"
colnames(IndexMx) <- names(mysets)
IndexMx <- data.frame(ProbesetID = annot$ProbesetID, IndexMx)

resAll <- merge(resAll, IndexMx, by = "ProbesetID", sort = FALSE)

write.table(resAll, file = "Comp1_DE_results_AllPlus.xls", quote = FALSE, sep = "\t", row.names = FALSE)


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups))

design <- model.matrix(~ 0 + Treatment, data=treatments)
rownames(design) <- targets$labels
design

contrasts <- cbind(CtrlCD4 = c(-1, 0, 0, 0, 1), CtrlCD4CD8 = c(0, -1, 0, 0, 1), CtrlCD8 = c(0, 0, -1, 0, 1), CtrlBM = c(0, 0, 0, -1, 1)) # treatment - control
contrasts




### run CAMERA

gsea <- list()

coef <- "CtrlCD4"
gsea.tmp <- gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=FALSE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("NGenes","Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), NGenes = gsea[[coef]][,1], gsea[[coef]][,-1])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)



### using information from eBayes fitting: fit2

pdf(paste0("PLOTS/GS_barcodeplot_",coef,".pdf"))

topgs <- 1
gsn <-rownames(gsea[[coef]])[1:topgs] 
gss <- gsea.tmp[gsn, , drop = FALSE]

for(i in 1:length(topgs)){
  
  barcodeplot(statistics = as.numeric((fit2$t[,coef])), index = Index[[gsn[i]]], index2 = NULL, gene.weights = as.numeric((fit2$coefficients[, coef]))[Index[[gsn[i]]]], weights.label = "logFC", labels = c("Up","Down"), quantiles = c(-1,1), col.bars = NULL, worm = TRUE, span.worm=0.45, main = paste0(gsn[i], "\n", gss[i, "Direction"], ", FDR = ", sprintf("%.02e",gss[i, "FDR"])))
  
  barcodeplot(statistics = as.numeric((fit2$p.value[, coef])), index = Index[[gsn[i]]], index2 = NULL, gene.weights = as.numeric((fit2$coefficients[, coef]))[Index[[gsn[i]]]], weights.label = "logFC", labels = c("Not significant","Significant"), quantiles = c(0.05,1), col.bars = NULL, worm = TRUE, span.worm=0.45, main = paste0(gsn[i], "\n", gss[i, "Direction"], ", FDR = ", sprintf("%.02e",gss[i, "FDR"])))
  
  
}

dev.off()




coef <- "CtrlCD4CD8"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)


coef <- "CtrlCD8"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)



coef <- "CtrlBM"
gsea[[coef]] <- camera(y = eset_main, index=Index, design=design, contrast=contrasts[,coef], trend.var=TRUE)
head(gsea[[coef]], 10)
table(gsea[[coef]]$FDR < 0.05)
gsea[[coef]] <- gsea[[coef]][, c("Direction", "PValue", "FDR")]
colnames(gsea[[coef]]) <- paste0(coef, "_", colnames(gsea[[coef]]))
gsea[[coef]] <- data.frame(GeneSet = rownames(gsea[[coef]]), gsea[[coef]])
# write.table(gsea[[coef]], paste("Comp1_GSEA_",coef ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)





### merge all results into one table
gseaAll <- merge(gsea[["CtrlCD4"]], gsea[["CtrlCD4CD8"]], by = "GeneSet", all = TRUE)
gseaAll <- merge(gseaAll, gsea[["CtrlCD8"]], by = "GeneSet", all = TRUE)
gseaAll <- merge(gseaAll, gsea[["CtrlBM"]], by = "GeneSet", all = TRUE)
write.table(gseaAll, paste("Comp1_GSEA_Hallmark_All.xls", sep=""), sep="\t", row.names=F, quote = FALSE)





# http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm




###########################################################################
### Clustering for all genes based on DE results (-1, 0, 1)
###########################################################################

targets <- targets_org
eset_main <- eset_main_org

### keep only leukemia and control CD4+, CD4+CD8+ and CD8+ samples
samples2keep <- targets_org$ExperimentShort != "afterTreatment" & targets_org$labels != "control_HeLa" & targets_org$labels != "control_wholeBoneMarrow"

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset_main <- eset_main[ ,ord]


expr <- exprs(eset_main)
### normalize expression per gene
exprNorm <- t(scale(t(expr), center = TRUE, scale = TRUE))




####### load the DE results

## does not work 
# resAll <- read.table("Comp1_DE_results_All.xls", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

library(limma)


allLines <- readLines("Comp1_DE_results_All.xls", n = -1)[-1]
resAll <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(resAll) <- strsplit2(readLines("Comp1_DE_results_All.xls", n = 1), "\t")

resAll <- resAll[, !grepl(pattern = "CEL", colnames(resAll))]

resAllSort <- resAll[order(resAll$CtrlCD4_PValues, resAll$CtrlCD4CD8_PValues, resAll$CtrlCD8_PValues, decreasing = FALSE), ]

resAllSort$clusters <- apply(resAllSort[, c("CtrlCD4_Results", "CtrlCD4CD8_Results", "CtrlCD8_Results")], MARGIN = 1, paste, collapse = ",")

resAllSort <- resAllSort[resAllSort$clusters != "0,0,0", ]


library(gtools)

clusters <- apply(permutations(n=3, r=3, v = c(-1, 0, 1), repeats.allowed=TRUE), MARGIN = 1, paste, collapse = ",")
clusters <- clusters[clusters != "0,0,0"]

resAllSort$clusters <- factor(resAllSort$clusters, levels = clusters)

resAllSort <- resAllSort[order(resAllSort$clusters), ]

# unique(resAllSort$clusters)

### number of genes in clusters
table(resAllSort$clusters)




##### Create a heat map with all the clusters

intrProbes <- as.character(resAllSort$ProbesetID)

# dataHeat <- expr[intrProbes, ]
dataHeat <- exprNorm[intrProbes, ]

annotation_col <- targets[, "groups", drop = FALSE]
rownames(annotation_col) <- colnames(dataHeat)

cols <- unique(targets$colors)
names(cols) <- unique(targets$group)

annotation_colors = list(groups = cols)

labels_row <- strsplit2(resAllSort$GeneSymbol, " /// ")[, 1]

library(pheatmap)


pdf("PLOTS/heatmap_clusters.pdf", width = 7, height = 10)

pheatmap(dataHeat, color = colorRamps::matlab.like(100), cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = rep("", nrow(dataHeat)), annotation_legend = FALSE, fontsize_row = 8, gaps_col = cumsum(table(targets$groups)), gaps_row = cumsum(table(resAllSort$clusters)),breaks =  seq(from = -4, to = 4, length.out = 101), legend_breaks = seq(from = -4, to = 4, by = 2))

dev.off()



write.table(resAllSort, file = "Comp1_DEclusters.xls", quote = FALSE, sep = "\t", row.names = FALSE)




###########################################################################
#### GO analysis per cluster 
###########################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")

library(topGO)
library(Rgraphviz)


affyLib <- "mogene20sttranscriptcluster.db"

### Function used to create new topGOdata object
fun.gene.sel <- function(geneList) {
  return(geneList <- ifelse(geneList==0, FALSE, TRUE))
}

### keep the clusters with at least 50 genes
cls <- levels(resAllSort$clusters)[table(resAllSort$clusters) > 50]


allResList <- list()

for(cl in cls){
  # cl <- cls[1]
  
  
  geneList <- rep(0, nrow(resAll))
  names(geneList) <- resAll$ProbesetID
  geneList[resAllSort[resAllSort$clusters == cl, "ProbesetID"]] <- 1
  table(geneList)
  
  
  for(go in c("BP","MF","CC")){
    # go = "BP"
    
    cat("Cluster:", cl, "go:", go, "\n")
    
    sampleGOdata <- new("topGOdata", description = paste0("Simple session for ", cl), ontology = go, allGenes = geneList, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
    
    #     print(sampleGOdata)
    
    result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    
    pValues <- score(result)
    topNodes <- length(pValues)
    
    allRes <- GenTable(sampleGOdata, elimFisher = result, orderBy = "elimFisher", topNodes = topNodes)      
    colnames(allRes)[6] <- "PValues" 
    allRes$GO <- go
    
    
    #     pdf(paste("PLOTS/GO_",cl, "_" ,go, ".pdf", sep=""))
    #     showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 5, useInfo = 'all')
    #     dev.off()
    
    allRes$AdjPValues <- p.adjust(allRes$PValues, method = "BH")
    
    #     cat("#########################################################################################", fill = TRUE)
    #     print(table(allRes$AdjPValues < 0.05))
    #     print(head(allRes, 20))
    #     cat("#########################################################################################", fill = TRUE)
    
    
    # write.table(allRes, paste("Comp1_GO_Fisher_elim_",cl, "_", go ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)
    
    allResList[[paste0(cl, "_", go)]] <- allRes
    
    
  }
  
}


save(allResList, file = "Comp1_GO_Clusters_Fisher_elim.rdata")

#### save results 

for(go in c("BP","MF","CC")){
  
  cl <- cls[1]
  allR <- allResList[[paste0(cl, "_", go)]]
  allAll <- allR[, c("GO.ID", "GO", "Term", "Annotated")]  
  
  for(cl in cls){
    # cl = cls[1]
    
    allR <- allResList[[paste0(cl, "_", go)]][, c("GO.ID", "Significant", "Expected", "PValues", "AdjPValues")] 
    ### add cluster names to columns
    colnames(allR) <- paste0(c("", rep(paste0("CL(",cl, ")_"), 4)), colnames(allR))
    
    ### merge all results into one table
    allAll <- merge(allAll, allR, by = "GO.ID", sort = FALSE)
    
  }
  
  
  write.table(allAll, paste0("Comp1_GO_Clusters_Fisher_elim_", go ,".xls"), sep="\t", row.names=F, quote = FALSE)
  
}




###########################################################################
#### GO analysis per control - up or down regulation 
###########################################################################

library(topGO)
library(Rgraphviz)


affyLib <- "mogene20sttranscriptcluster.db"

### Function used to create new topGOdata object
fun.gene.sel <- function(geneList) {
  return(geneList <- ifelse(geneList==0, FALSE, TRUE))
}


cls <- rep(c(-1, 1), times = 4)
names(cls) <- rep(c("CtrlCD4", "CtrlCD4CD8", "CtrlCD8", "CtrlBM"), each = 2)

allResList <- list()

for(cl in 1:length(cls)){
  # cl <- cls[1]
  
  cl <- cls[cl]
  
  print(cl)
  geneList <- rep(0, nrow(resAll))
  names(geneList) <- resAll$ProbesetID
  
  geneList[resAll[resAll[, paste0(names(cl), "_Results")] == cl, "ProbesetID"]] <- 1
  
  table(geneList)
  
  cl <- paste0(names(cl),".", cl)
  
  for(go in c("BP","MF","CC")){
    # go = "BP"
    
    cat("Cluster:", cl, "go:", go, "\n")
    
    sampleGOdata <- new("topGOdata", description = paste0("Simple session for ", cl), ontology = go, allGenes = geneList, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
    
    #     print(sampleGOdata)
    
    result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    
    pValues <- score(result)
    topNodes <- length(pValues)
    
    allRes <- GenTable(sampleGOdata, elimFisher = result, orderBy = "elimFisher", topNodes = topNodes)      
    colnames(allRes)[6] <- "PValues" 
    allRes$GO <- go
    
    
    #     pdf(paste("PLOTS/GO_",cl, "_" ,go, ".pdf", sep=""))
    #     showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 5, useInfo = 'all')
    #     dev.off()
    
    allRes$AdjPValues <- p.adjust(allRes$PValues, method = "BH")
    
    #     cat("#########################################################################################", fill = TRUE)
    #     print(table(allRes$AdjPValues < 0.05))
    #     print(head(allRes, 20))
    #     cat("#########################################################################################", fill = TRUE)
    
    
    # write.table(allRes, paste("Comp1_GO_Fisher_elim_",cl, "_", go ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)
    
    allResList[[paste0(cl, "_", go)]] <- allRes
    
    
  }
  
}


save(allResList, file = "Comp1_GO_UpDown_Fisher_elim.rdata")


#### save results 


cls <- paste0(names(cls), "." ,cls)

for(go in c("BP","MF","CC")){
  
  cl <- cls[1]
  
  allR <- allResList[[paste0(cl, "_", go)]]
  allAll <- allR[, c("GO.ID", "GO", "Term", "Annotated")]  
  
  for(cl in cls){
    # cl = cls[1]
    
    allR <- allResList[[paste0(cl, "_", go)]][, c("GO.ID", "Significant", "Expected", "PValues", "AdjPValues")] 
    ### add cluster names to columns
    colnames(allR) <- paste0(c("", rep(paste0("CL(",cl, ")_"), 4)), colnames(allR))
    
    ### merge all results into one table
    allAll <- merge(allAll, allR, by = "GO.ID", sort = FALSE)
    
  }
  
  
  write.table(allAll, paste0("Comp1_GO_UpDown_Fisher_elim_", go ,".xls"), sep="\t", row.names=F, quote = FALSE)
  
}




###########################################################################
#### GO analysis per control 
###########################################################################

# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")

library(topGO)
library(Rgraphviz)


affyLib <- "mogene20sttranscriptcluster.db"

### Function used to create new topGOdata object
fun.gene.sel <- function(geneList) {
  return(geneList <- ifelse(geneList==0, FALSE, TRUE))
}

coefs <- c("CtrlCD4", "CtrlCD4CD8", "CtrlCD8", "CtrlBM")

allResList <- list()

for(coef in coefs){
  # coef <- coefs[1]
  
  tt <- topTable(fit2, coef=coef, n=Inf)
  geneList <- rep(0, nrow(tt))
  names(geneList) <- rownames(tt)
  geneList[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1] <- 1
  print(table(geneList))
  
  for(go in c("BP","MF","CC")){
    # go = "BP"
    print(coef)
    print(go)
    sampleGOdata <- new("topGOdata", description = paste0("Simple session for ", coef), ontology = go, allGenes = geneList, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
    
    print(sampleGOdata)
    
    result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    
    pValues <- score(result)
    topNodes <- length(pValues)
    
    allRes <- GenTable(sampleGOdata, elimFisher = result, orderBy = "elimFisher", topNodes = topNodes)      
    colnames(allRes)[6] <- "PValues" 
    
    pdf(paste("PLOTS/GO_",coef, "_" ,go, ".pdf", sep=""))
    showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 5, useInfo = 'all')
    dev.off()
    
    allRes$AdjPValues <- p.adjust(allRes$PValues, method = "BH")
    
    print(table(allRes$AdjPValues < 0.05))
    
    cat("#########################################################################################", fill = TRUE)
    print(head(allRes, 20))
    cat("#########################################################################################", fill = TRUE)
    
    
    # write.table(allRes, paste("Comp1_GO_Fisher_elim_",coef, "_", go ,".xls", sep=""), sep="\t", row.names=F, quote = FALSE)
    
    allResList[[paste0(go, "_", coef)]] <- allRes
    
    
  }
  
}


#### save results 

for(go in c("BP","MF","CC")){
  
  coef <- "CtrlCD4"
  allR <- allResList[[paste0(go, "_", coef)]]
  colnames(allR) <- paste0(c(rep("", 3), rep(paste0(coef, "_"), 4)), colnames(allR))
  allResList[[paste0(go, "_", coef)]] <- allR
  
  coef <- "CtrlCD4CD8"
  allR <- allResList[[paste0(go, "_", coef)]][, -c(2, 3)]
  colnames(allR) <- paste0(c("", rep(paste0(coef, "_"), 4)), colnames(allR))
  allResList[[paste0(go, "_", coef)]] <- allR
  
  
  coef <- "CtrlCD8"
  allR <- allResList[[paste0(go, "_", coef)]][, -c(2, 3)]
  colnames(allR) <- paste0(c("", rep(paste0(coef, "_"), 4)), colnames(allR))
  allResList[[paste0(go, "_", coef)]] <- allR
  
  
  coef <- "CtrlBM"
  allR <- allResList[[paste0(go, "_", coef)]][, -c(2, 3)]
  colnames(allR) <- paste0(c("", rep(paste0(coef, "_"), 4)), colnames(allR))
  allResList[[paste0(go, "_", coef)]] <- allR
  
  
  ### merge all results into one table
  allAll <- merge(allResList[[paste0(go, "_", "CtrlCD4")]], allResList[[paste0(go, "_", "CtrlCD4CD8")]], by = "GO.ID", all = TRUE)  
  allAll <- merge(allAll, allResList[[paste0(go, "_", "CtrlCD8")]], by = "GO.ID", all = TRUE)
  allAll <- merge(allAll, allResList[[paste0(go, "_", "CtrlBM")]], by = "GO.ID", all = TRUE)
  write.table(allAll, paste0("Comp1_GO_Fisher_elim_", go ,".xls"), sep="\t", row.names=F, quote = FALSE)
  
}






###########################################################################
#### Comparison 2: ALL pre VS ALL after treatment
###########################################################################


targets <- targets_org
eset_main <- eset_main_org

### keep only leukemia and afterTreatment samples
samples2keep <- grepl("leukemia|afterTreatment", targets$labels)

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset_main <- eset_main[ ,ord]

# all(sampleNames(eset_main) == strsplit2(targets$FileName, "//")[,2])


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups))
treatments$Treatment <- relevel(treatments$Treatment, ref = "leukemia")
treatments


design <- model.matrix(~Treatment, data=treatments)
rownames(design) <- targets$labels
design



fit <- lmFit(eset_main, design)

fit2 <- eBayes(fit[, "TreatmentafterTreatment"], trend = TRUE)


pdf(paste0(path_plots, "Comp2_plotSA_trend.pdf"))
plotSA(fit2)
dev.off()


## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)


colours <- unique(targets[targets$groups == "afterTreatment", "colors"])

pdf(paste0(path_plots, "Comp2_vennDiagram.pdf"))
vennDiagram(results,include=c("up", "down"), circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="both", circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="up", circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="down", circle.col=colours, counts.col=c("gold", "darkblue"))
dev.off()




table <- topTable(fit2, coef = 1, n = Inf)


### save all results with nice order
resCoeff <- fit2$coefficients
resT <- fit2$t
resPValue <- fit2$p.value
resPValueAdj <- apply(fit2$p.value, 2, p.adjust, method = "BH")
resRes <- results[, 1]

resDE <- data.frame(resCoeff, resT, resPValue, resPValueAdj, resRes)
colnames(resDE) <- paste0("afterTreatment_all_", c("coeffs", "t", "PValues", "AdjPValues", "Results"))

resGenes <- fit2$genes

resExpr <- round(exprs(eset_main_org), 2)
colnames(resExpr) <- paste0(targets_org$labels, "_", colnames(resExpr))
resExpr <- resExpr[, order(colnames(resExpr))]

resAll <- cbind(resGenes, resDE, resExpr)

write.table(resAll, file = paste0(path_results, "Comp2_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)


### plot MA
pdf(paste0(path_plots, "Comp2_plotMA.pdf"))

limma::plotMA(fit2, coef = 1, status = results, values = c(-1, 0, 1), col = c("red", "black", "green"), cex = c(0.7, 0.3, 0.7))
abline(0,0,col="blue")

dev.off()




### volcano plots
library(ggplot2)


table <- topTable(fit2, coef = 1, n=Inf)
table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
gg2 <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + theme_bw() + theme(legend.position = "none") +  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("after Treatment")

pdf(paste0(path_plots, "Comp2_volcanoplot.pdf"))
print(gg2)
dev.off()



### histograms of p-values and adjusted p-values
colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp2_hist_pvs.pdf"))

table <- topTable(fit2, coef = 1, n=Inf)
hist(table$P.Value, breaks = 100, main = "afterTreatment", xlab = "P-values", col = colours)

dev.off()




###########################################################################
#### Comparison 3a:  pre VS after treatment with matched samples pooled
###########################################################################

library(oligo)
library(pd.mogene.2.0.st)
library(limma)

load(paste0(path_results, "eset_main_org.Rdata"))
load(paste0(path_results, "targets_org.Rdata"))

targets <- targets_org
eset_main <- eset_main_org

tt <- table(targets$CellTypeShort, targets$groups)

tt <- data.frame(tt, stringsAsFactors = FALSE)

cell_types <- as.character(tt[tt$Var2 == "afterTreatment" & tt$Freq > 0, "Var1"])

### keep only leukemia and afterTreatment samples that have matched cell type 
samples2keep <- grepl("leukemia|afterTreatment", targets$labels) & targets$CellTypeShort %in% cell_types

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset_main <- eset_main[ ,ord]

# all(sampleNames(eset_main) == strsplit2(targets$FileName, "//")[,2])


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups), CellType = targets$CellTypeShort)
treatments$Treatment <- relevel(treatments$Treatment, ref = "leukemia")
treatments


design <- model.matrix(~ Treatment, data = treatments)
rownames(design) <- targets$labels
design



fit <- lmFit(eset_main, design)

fit2 <- eBayes(fit[, "TreatmentafterTreatment"], trend = TRUE)


pdf(paste0(path_plots, "Comp3a_plotSA_trend.pdf"))
plotSA(fit2)
dev.off()


## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)


colours <- unique(targets[targets$groups == "afterTreatment", "colors"])

pdf(paste0(path_plots, "Comp3a_vennDiagram.pdf"))
vennDiagram(results,include=c("up", "down"), circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="both", circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="up", circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="down", circle.col=colours, counts.col=c("gold", "darkblue"))
dev.off()




table <- topTable(fit2, coef = 1, n = Inf)


### save all results with nice order
resCoeff <- fit2$coefficients
resT <- fit2$t
resPValue <- fit2$p.value
resPValueAdj <- apply(fit2$p.value, 2, p.adjust, method = "BH")
resRes <- results[, 1]

resDE <- data.frame(resCoeff, resT, resPValue, resPValueAdj, resRes)
colnames(resDE) <- paste0("afterTreatment_matched_pooled_", c("coeffs", "t", "PValues", "AdjPValues", "Results"))

resGenes <- fit2$genes

resExpr <- round(exprs(eset_main_org), 2)
colnames(resExpr) <- paste0(targets_org$labels, "_", colnames(resExpr))
resExpr <- resExpr[, order(colnames(resExpr))]

resAll <- cbind(resGenes, resDE, resExpr)

write.table(resAll, file = paste0(path_results, "Comp3a_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




### plot MA
pdf(paste0(path_plots, "Comp3a_plotMA.pdf"))

limma::plotMA(fit2, coef = 1, status = results, values = c(-1, 0, 1), col = c("red", "black", "green"), cex = c(0.7, 0.3, 0.7))
abline(0,0,col="blue")

dev.off()




### volcano plots
library(ggplot2)


table <- topTable(fit2, coef = 1, n=Inf)
table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
gg2 <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + theme_bw() + theme(legend.position = "none") +  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("after Treatment")

pdf(paste0(path_plots, "Comp3a_volcanoplot.pdf"))
print(gg2)
dev.off()



### histograms of p-values and adjusted p-values
colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp3a_hist_pvs.pdf"))

table <- topTable(fit2, coef = 1, n=Inf)
hist(table$P.Value, breaks = 100, main = "afterTreatment", xlab = "P-values", col = colours)

dev.off()




### plot expression of top sign. genes/probesets
library(ggplot2)
library(reshape2)

expr <- exprs(eset_main)
topn <- 20

rownames(targets) <- strsplit2(targets$FileName, split = "//")[, 2] 

tt <- topTable(fit2, coef = 1, n=Inf, p.value=0.05, lfc=1)

### in the report display only first gene symbol
GeneSymbol <- strsplit2(head(tt[,"GeneSymbol"], topn), " /// ")[,1]
GeneTitle <- paste0(substr(strsplit2(head(tt[,"GeneTitle"], topn), " /// ")[,1], 1, 30))

# print(data.frame(GeneSymbol = GeneSymbol, GeneTitle = GeneTitle , head(tt[, c("logFC", "AveExpr", "P.Value", "adj.P.Val")], topn)))

topp <- rownames(tt)[1:topn]

df <- data.frame(Gene = topp, expr[topp,])
df.m <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")
### keep order of genes as in tt
df.m$Gene <- factor(df.m$Gene, levels = topp)
### add Entrez ID to the facet labels
lab.fct <- paste0(topp, "\n", strsplit2(tt[topp, "GeneSymbol"], " /// ")[,1])
levels(df.m$Gene) <- lab.fct
df.m$groups <- targets[df.m$Sample ,"groups"]

fill_colors <- unique(targets[, c("groups", "colors")])
fill_colors <- fill_colors[order(fill_colors$groups), "colors"]

ggp <- ggplot(df.m, aes(x = Sample, y = Expression, fill = groups)) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 10), plot.title = element_text(size = 16), strip.text.x = element_text(size = 10)) +
  scale_x_discrete(labels=targets$CellTypeShort) +
  labs(y = "Log2 expression") +
  geom_bar(stat = "identity") +
  facet_wrap(~ Gene, scales="free_y", ncol=4) +
  scale_fill_manual(values = fill_colors)

pdf(paste0(path_plots, "Comp3a_topExpressionBarPlot.pdf"), 11, 11)
print(ggp)    
dev.off()



###########################################################################
#### Comparison 3b:  pre VS after treatment with matched samples + cell type
###########################################################################

library(oligo)
library(pd.mogene.2.0.st)
library(limma)

load(paste0(path_results, "eset_main_org.Rdata"))
load(paste0(path_results, "targets_org.Rdata"))

targets <- targets_org
eset_main <- eset_main_org

tt <- table(targets$CellTypeShort, targets$groups)

tt <- data.frame(tt, stringsAsFactors = FALSE)

cell_types <- as.character(tt[tt$Var2 == "afterTreatment" & tt$Freq > 0, "Var1"])

### keep only leukemia and afterTreatment samples that have matched cell type 
samples2keep <- grepl("leukemia|afterTreatment", targets$labels) & targets$CellTypeShort %in% cell_types

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset_main <- eset_main[ ,ord]

# all(sampleNames(eset_main) == strsplit2(targets$FileName, "//")[,2])


#### design & analysis

treatments <- data.frame(Treatment = as.character(targets$groups), CellType = targets$CellTypeShort)
treatments$Treatment <- relevel(treatments$Treatment, ref = "leukemia")
treatments


design <- model.matrix(~ 0 + CellType + Treatment, data = treatments)
rownames(design) <- targets$labels
design



fit <- lmFit(eset_main, design)

fit2 <- eBayes(fit[, "TreatmentafterTreatment"], trend = TRUE)


pdf(paste0(path_plots, "Comp3b_plotSA_trend.pdf"))
plotSA(fit2)
dev.off()


## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)


colours <- unique(targets[targets$groups == "afterTreatment", "colors"])

pdf(paste0(path_plots, "Comp3b_vennDiagram.pdf"))
vennDiagram(results,include=c("up", "down"), circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="both", circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="up", circle.col=colours, counts.col=c("gold", "darkblue"))
# vennDiagram(results,include="down", circle.col=colours, counts.col=c("gold", "darkblue"))
dev.off()




table <- topTable(fit2, coef = 1, n = Inf)


### save all results with nice order
resCoeff <- fit2$coefficients
resT <- fit2$t
resPValue <- fit2$p.value
resPValueAdj <- apply(fit2$p.value, 2, p.adjust, method = "BH")
resRes <- results[, 1]

resDE <- data.frame(resCoeff, resT, resPValue, resPValueAdj, resRes)
colnames(resDE) <- paste0("afterTreatment_matched_paired_", c("coeffs", "t", "PValues", "AdjPValues", "Results"))

resGenes <- fit2$genes

resExpr <- round(exprs(eset_main_org), 2)
colnames(resExpr) <- paste0(targets_org$labels, "_", colnames(resExpr))
resExpr <- resExpr[, order(colnames(resExpr))]

resAll <- cbind(resGenes, resDE, resExpr)

write.table(resAll, file = paste0(path_results, "Comp3b_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




### plot MA
pdf(paste0(path_plots, "Comp3b_plotMA.pdf"))

limma::plotMA(fit2, coef = 1, status = results, values = c(-1, 0, 1), col = c("red", "black", "green"), cex = c(0.7, 0.3, 0.7))
abline(0,0,col="blue")

dev.off()




### volcano plots
library(ggplot2)


table <- topTable(fit2, coef = 1, n=Inf)
table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
gg2 <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + theme_bw() + theme(legend.position = "none") +  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("after Treatment")

pdf(paste0(path_plots, "Comp3b_volcanoplot.pdf"))
print(gg2)
dev.off()



### histograms of p-values and adjusted p-values
colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp3b_hist_pvs.pdf"))

table <- topTable(fit2, coef = 1, n=Inf)
hist(table$P.Value, breaks = 100, main = "afterTreatment", xlab = "P-values", col = colours)

dev.off()




### plot expression of top sign. genes/probesets
library(ggplot2)
library(reshape2)

expr <- exprs(eset_main)
topn <- 20

rownames(targets) <- strsplit2(targets$FileName, split = "//")[, 2] 

tt <- topTable(fit2, coef = 1, n=Inf, p.value=0.05, lfc=1)

### in the report display only first gene symbol
GeneSymbol <- strsplit2(head(tt[,"GeneSymbol"], topn), " /// ")[,1]
GeneTitle <- paste0(substr(strsplit2(head(tt[,"GeneTitle"], topn), " /// ")[,1], 1, 30))

# print(data.frame(GeneSymbol = GeneSymbol, GeneTitle = GeneTitle , head(tt[, c("logFC", "AveExpr", "P.Value", "adj.P.Val")], topn)))

topp <- rownames(tt)[1:topn]

df <- data.frame(Gene = topp, expr[topp,])
df.m <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")
### keep order of genes as in tt
df.m$Gene <- factor(df.m$Gene, levels = topp)
### add Entrez ID to the facet labels
lab.fct <- paste0(topp, "\n", strsplit2(tt[topp, "GeneSymbol"], " /// ")[,1])
levels(df.m$Gene) <- lab.fct
df.m$groups <- targets[df.m$Sample ,"groups"]

fill_colors <- unique(targets[, c("groups", "colors")])
fill_colors <- fill_colors[order(fill_colors$groups), "colors"]

ggp <- ggplot(df.m, aes(x = Sample, y = Expression, fill = groups)) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 10), plot.title = element_text(size = 16), strip.text.x = element_text(size = 10)) +
  scale_x_discrete(labels=targets$CellTypeShort) +
  labs(y = "Log2 expression") +
  geom_bar(stat = "identity") +
  facet_wrap(~ Gene, scales="free_y", ncol=4) +
  scale_fill_manual(values = fill_colors)

pdf(paste0(path_plots, "Comp3b_topExpressionBarPlot.pdf"), 11, 11)
print(ggp)    
dev.off()









###########################################################################
#### Merge all results
###########################################################################


res_files <- c("Comp1_DE_results_All.xls", "Comp2_DE_results_All.xls", "Comp3a_DE_results_All.xls", "Comp3b_DE_results_All.xls")


res_all <- lapply(1:length(res_files), function(ff){
  # ff = 1
  
  allLines <- readLines(paste0(path_results, res_files[ff]), n = -1)[-1]
  resComp <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
  colnames(resComp) <- strsplit2(readLines(paste0(path_results, res_files[ff]), n = 1), "\t")
  
  
  if(ff == 1){
    return(resComp[, !grepl(pattern = "CEL", x = colnames(resComp))])
  }else if(ff == length(res_files)){
    return(resComp[, !grepl(pattern = "Gene", x = colnames(resComp))])
  }else{
    return(resComp[, !grepl(pattern = "Gene", x = colnames(resComp)) & !grepl(pattern = "CEL", x = colnames(resComp))])
  }
  
})

lapply(res_all, colnames)


res_all <- Reduce(function(...) merge(..., by = "ProbesetID", all=TRUE, sort = FALSE), res_all)

colnames(res_all)

write.table(res_all, file = paste0(path_results, "CompALL_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




results <- res_all[, grepl("afterTreatment.*Results", colnames(res_all))]

colnames(results) <- gsub("_Results", "", gsub(pattern = "afterTreatment_", "", colnames(results)))


pdf(paste0(path_plots, "CompAll_afterTreatment_vennDiagram.pdf"))
vennDiagram(results, include = c("up", "down"), counts.col = c("gold", "darkblue"))
dev.off()































