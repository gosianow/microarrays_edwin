###########################################################################
# Created 23 Nov 2015
# BioC 3.1

# DE analysis of Affymetrix Mouse Gene 2.0 ST arrays (pd.mogene.2.0.st)
# Additional replicates 

# Pre versus after treatment analysis 

# Modified 10 Dec 2015
# GO analysis 

###########################################################################

setwd("/home/Shared/data/array/microarrays_edwin")

path_plots <- "Analysis_Mouse_Rebuttal_v2/Plots/"
path_results <- "Analysis_Mouse_Rebuttal_v2/"

dir.create(path_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(path_results, showWarnings = FALSE, recursive = TRUE)

###########################################################################

library(limma)
library(plyr)
library(oligo)
library(pd.mogene.2.0.st)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(topGO)
library(Rgraphviz)
library(mogene20sttranscriptcluster.db)


###########################################################################
# create targets table with information about samples from Micro_array_sample_list.txt
###########################################################################


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
targets_batch1$batch <- "batch1"

targets_batch2 <- read.table(file.path("metadata", "Micro_array_sample_list_rebuttal.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)
targets_batch2$batch <- "batch2"

colors <- unique(targets_batch1[, c("groups", "colors")])
targets_batch2$colors <- colors$colors[match(targets_batch2$groups, colors$groups)] 

targets <- rbind.fill(targets_batch1, targets_batch2)

targets$TypeOfReplicate[is.na(targets$TypeOfReplicate)] <- "biological"

targets$TypeOfReplicate[targets$CellTypeShort == "B2M3" & targets$groups == "afterTreatment"] <- "technical"


write.table(targets, file = file.path("metadata", "targets_all.xls"), quote = FALSE, sep = "\t", row.names = FALSE)


###########################################################################
# read in all targets
###########################################################################


targets_org <- targets <- read.table(file.path("metadata", "targets_all.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)

head(targets)

###########################################################################
#### import cel files
###########################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("pd.mogene.2.0.st")


ff <- as.character(targets$FileName)

x <- oligo::read.celfiles(filenames = ff) ## GeneFeatureSet

pdf(paste0(path_plots, "PreProc_boxplot.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
boxplot(x, las = 2, col = targets$colors, names = targets$labels, las=2)
dev.off()

pdf(paste0(path_plots, "PreProc_hist.pdf"))
par(mar = c(5, 4, 4, 2) + 0.1)
hist(x, col = targets$colors, lwd = 2)
legend("topright", legend = targets$labels, col =  targets$colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()



###########################################################################
### PLM normalization; create images of chips, NUSE and RLE plots
###########################################################################


fitplm <- oligo::fitProbeLevelModel(x)


pdf(paste0(path_plots,"PreProc_NUSE_fitplm.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
oligo::NUSE(fitplm, col = targets$colors, names = targets$labels, las=2)
dev.off()


pdf(paste0(path_plots, "PreProc_RLE_fitplm.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
oligo::RLE(fitplm, col = targets$colors, names = targets$labels, las=2)
dev.off()



###########################################################################
### Normalization with RMA
###########################################################################


eset_org <- eset <- oligo::rma(x) ## Is the expression in log2 scale? ## ExpressionSet


pdf(paste0(path_plots, "PreProc_boxplot_norm.pdf"))
par(mar = c(12, 4, 4, 2) + 0.1) # c(bottom, left, top, right), default = c(5, 4, 4, 2) + 0.1
boxplot(eset, las = 2, col = targets$colors, names = targets$labels)
dev.off()

pdf(paste0(path_plots, "PreProc_hist_norm.pdf"))
par(mar = c(5, 4, 4, 2) + 0.1)
hist(eset, col = targets$colors, lwd = 2)
legend("topright", legend = targets$labels, col =  targets$colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()



###########################################################################
### MDS plots
###########################################################################

########## All samples

eset <- eset_org
targets <- targets_org

labels <- targets$groups


pdf(paste0(path_plots, "MDS_all.pdf"), width = 5, height = 5)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
dev.off()


legend <- unique(targets[, c("groups", "colors")])
rownames(legend) <- legend$groups

rx <- range(mds$x)
ry <- range(mds$y)

rplot <- max(rx[2] - rx[1], ry[2] - ry[1])/2

xmin <- mean(rx) - rplot
xmax <- mean(rx) + rplot
ymin <- mean(ry) - rplot
ymax <- mean(ry) + rplot



# pdf(paste0(path_plots, "MDS_all_points.pdf"), width = 5, height = 5)
# plot(mds$x, mds$y, pch = as.numeric(factor(targets$batch)), col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
# # text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.3, cex = 0.3)
# legend("bottomleft", legend = c("batch1", "batch2", legend$groups), pch = c(1, 2, rep(16, nrow(legend))), col = c(1, 1, legend$colors), cex = 0.8, bty = "n")
# dev.off()



mdsdf <- data.frame(x = mds$x, y = mds$y, groups = factor(targets$groups, levels = legend$groups), batch = targets$batch, TypeOfReplicate = targets$TypeOfReplicate, CellTypeShort = targets$CellTypeShort)


ggp <- ggplot(mdsdf, aes(x = x, y = y, colour = groups, shape = batch)) +
  theme_bw() +
  geom_point(size = 4) +
  geom_point(data = mdsdf, aes(x = x, y = y, fill = groups, alpha = TypeOfReplicate), size = 4) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "right", panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_shape_manual(values = c(21, 24)) +
  scale_colour_manual(values = legend$colors) +
  scale_fill_manual(values = legend$colors) +
  scale_alpha_manual(values = c(0.05, 0.3)) +
  coord_cartesian(xlim = c(xmin-0.1, xmax+0.1), ylim = c(ymin-0.1, ymax+0.1))


pdf(paste0(path_plots, "MDS_all_points_gg.pdf"), width = 8, height = 5)
print(ggp)
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
rownames(legend) <- legend$groups

rx <- range(mds$x)
ry <- range(mds$y)

rplot <- max(rx[2] - rx[1], ry[2] - ry[1])/2

xmin <- mean(rx) - rplot
xmax <- mean(rx) + rplot
ymin <- mean(ry) - rplot
ymax <- mean(ry) + rplot


# pdf(paste0(path_plots, "MDS_all_noHela_points.pdf"), width = 5, height = 5)
# plot(mds$x, mds$y, pch = as.numeric(factor(targets$batch)), col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
# # text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.3, cex = 0.3)
# legend("bottomleft", legend = c("batch1", "batch2", legend$groups), pch = c(1, 2, rep(16, nrow(legend))), col = c(1, 1, legend$colors), cex = 0.8, bty = "n")
# dev.off()





mdsdf <- data.frame(x = mds$x, y = mds$y, groups = factor(targets$groups, levels = legend$groups), batch = targets$batch, TypeOfReplicate = targets$TypeOfReplicate, CellTypeShort = targets$CellTypeShort)

ggp <- ggplot(mdsdf, aes(x = x, y = y, colour = groups, shape = batch)) +
  theme_bw() +
  geom_point(size = 4) +
  geom_point(data = mdsdf, aes(x = x, y = y, fill = groups, alpha = TypeOfReplicate), size = 4) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "right", panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_shape_manual(values = c(21, 24)) +
  scale_colour_manual(values = legend$colors) +
  scale_fill_manual(values = legend$colors) +
  scale_alpha_manual(values = c(0.05, 0.3)) +
  coord_cartesian(xlim = c(xmin-0.1, xmax+0.1), ylim = c(ymin-0.1, ymax+0.1))


pdf(paste0(path_plots, "MDS_all_noHela_points_gg.pdf"), width = 8, height = 5)
print(ggp)
dev.off()




### zoom on leukemia and treatment samples

keep_samps <- targets$ExperimentShort != "control"

eset <- eset[, keep_samps]
targets <- targets[keep_samps, ]

mds$x <- mds$x[keep_samps]
mds$y <- mds$y[keep_samps]


legend <- unique(targets[, c("groups", "colors")])
rownames(legend) <- legend$groups

rx <- range(mds$x)
ry <- range(mds$y)

rplot <- max(rx[2] - rx[1], ry[2] - ry[1])/2

xmin <- mean(rx) - rplot
xmax <- mean(rx) + rplot
ymin <- mean(ry) - rplot
ymax <- mean(ry) + rplot



# pdf(paste0(path_plots, "MDS_all_noHela_points_zoom.pdf"), width = 5, height = 5)
# plot(mds$x, mds$y, pch = as.numeric(factor(targets$batch)), col = targets$colors, las = 1, cex.axis = 1, cex.lab = 1, xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2", cex = 1, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
# text(mds$x, mds$y, labels = targets$CellTypeShort, pos = 3, offset = 0.4, cex = 0.5, col = targets$colors)
# legend("bottomleft", legend = c("batch1", "batch2", legend$groups), pch = c(1, 2, rep(16, nrow(legend))), col = c(1, 1, legend$colors), cex = 0.8, bty = "n")
# dev.off()



mdsdf <- data.frame(x = mds$x, y = mds$y, groups = factor(targets$groups, levels = legend$groups), batch = targets$batch, TypeOfReplicate = targets$TypeOfReplicate, CellTypeShort = targets$CellTypeShort)


ggp <- ggplot(mdsdf, aes(x = x, y = y, colour = groups, shape = batch)) +
  theme_bw() +
  geom_point(size = 4) +
  geom_point(data = mdsdf, aes(x = x, y = y, fill = groups, alpha = TypeOfReplicate), size = 4) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  geom_text(aes(label = CellTypeShort), hjust = 1/2, vjust = -1, size = 3) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "right", panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_shape_manual(values = c(21, 24)) +
  scale_colour_manual(values = legend$colors) +
  scale_fill_manual(values = legend$colors) +
  scale_alpha_manual(values = c(0.05, 0.3)) +
  coord_cartesian(xlim = c(xmin-0.1, xmax+0.1), ylim = c(ymin-0.1, ymax+0.1))



pdf(paste0(path_plots, "MDS_all_noHela_points_zoom_gg.pdf"), width = 8, height = 5)
print(ggp)
dev.off()



########## Only leukemia (pre and post treatment) samples

keep_samps <- targets_org$groups %in% c("leukemia", "afterTreatment") 

eset <- eset_org[, keep_samps]
targets <- targets_org[keep_samps, ]

labels <- targets$groups


pdf(paste0(path_plots, "MDS_leukemia.pdf"), width = 5, height = 5)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)
dev.off()


legend <- unique(targets[, c("groups", "colors")])
rownames(legend) <- legend$groups

rx <- range(mds$x)
ry <- range(mds$y)

rplot <- max(rx[2] - rx[1], ry[2] - ry[1])/2

xmin <- mean(rx) - rplot
xmax <- mean(rx) + rplot
ymin <- mean(ry) - rplot
ymax <- mean(ry) + rplot



mdsdf <- data.frame(x = mds$x, y = mds$y, groups = factor(targets$groups, levels = legend$groups), batch = targets$batch, TypeOfReplicate = targets$TypeOfReplicate, CellTypeShort = targets$CellTypeShort)

ggp <- ggplot(mdsdf, aes(x = x, y = y, colour = groups, shape = batch)) +
  theme_bw() +
  geom_point(size = 4) +
  geom_point(data = mdsdf, aes(x = x, y = y, fill = groups, alpha = TypeOfReplicate), size = 4) +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  geom_text(aes(label = CellTypeShort), hjust = 1/2, vjust = -1, size = 3) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = "right", panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_shape_manual(values = c(21, 24)) +
  scale_colour_manual(values = legend$colors) +
  scale_fill_manual(values = legend$colors) +
  scale_alpha_manual(values = c(0.05, 0.3)) +
  coord_cartesian(xlim = c(xmin-0.1, xmax+0.1), ylim = c(ymin-0.1, ymax+0.1))


pdf(paste0(path_plots, "MDS_leukemia_points_gg.pdf"), width = 8, height = 5)
print(ggp)
dev.off()




###########################################################################
####### Do NOT keep the HeLa sample for the rest of the analysis
###########################################################################



keepSAMPS <- !targets_org$labels == "control_HeLa" 

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
####### Get annotation from formated files from Affy website (NetAffx Analysis Center)
### http://www.affymetrix.com/analysis/index.affx
###########################################################################

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

listMarts(host="www.ensembl.org")


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
listDatasets(mart)


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")

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

genes_merge <- plyr::ddply(genes, "affy_mogene_2_1_st_v1", plyr::summarize, GeneSymbol_Ensembl = paste0(external_gene_name, collapse = " /// "), GeneTitle_Ensembl = paste0(description, collapse = " /// "), EnsemblGeneID = paste0(ensembl_gene_id, collapse = " /// "))


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



all(annot_merge$ProbesetID == featureNames(eset_main)) ### is FALSE


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







bgTransInfo <- plyr::ddply(bgProbeInfo, "transcript_cluster_id", plyr::summarize,  gcCont=mean(gcCont))


bgdf <- data.frame(bgTransInfo, bgExpr)

bgdfm <- melt(bgdf, id.vars = c("transcript_cluster_id", "gcCont"), variable.name = "Samples", value.name = "Expression")
head(bgdfm)

bgdfm$gcCont <- factor(bgdfm$gcCont)

ggp.bg <- ggplot(data = bgdfm, aes(x = gcCont, y = Expression)) +
  geom_boxplot(colour = "lightcoral") +
  theme_bw()

pdf(paste0(path_plots, "PreProc_gc_boxplot.pdf"))
print(ggp.bg)
dev.off()


expr <- exprs(eset_main)

### Get the GC content for all the probe sets
transInfo <- plyr::ddply(probeInfo, "transcript_cluster_id", plyr::summarize,  gcCont=mean(gcCont))
rownames(transInfo) <- transInfo$transcript_cluster_id
transInfo <- transInfo[rownames(expr), ]
transInfo$gcCont <- round(transInfo$gcCont)
### see what is the average GC content for main probe sets
table(transInfo$gcCont)

df <- data.frame(transInfo, expr)

dfm <- melt(df, id.vars = c("transcript_cluster_id", "gcCont"), variable.name = "Samples", value.name = "Expression")
head(dfm)


dfm$Type <- "Main"
bgdfm$Type <- "BGP"

df.all <- rbind(dfm, bgdfm)
df.all$gcCont <- factor(df.all$gcCont, levels = 3:25)


ggp <- ggplot(data = df.all, aes(x = gcCont, y = Expression, fill = Type)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position="top")

pdf(paste0(path_plots, "PreProc_gc_boxplot_main_and_bgp.pdf"))
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
#### Comparison 3b:  pre VS after treatment with matched samples + cell type
### DO use the B2M3 technical replicate from batch 2 
###########################################################################

res_name <- "afterTreatment_matched_paired_withTechnical_"
comp_name <- "3b"

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
ord <- order(targets$groups, targets$CellTypeShort)
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


pdf(paste0(path_plots, "Comp", comp_name ,"_plotSA_trend.pdf"))
plotSA(fit2)
dev.off()


## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)


colours <- unique(targets[targets$groups == "afterTreatment", "colors"])

pdf(paste0(path_plots, "Comp", comp_name ,"_vennDiagram.pdf"))
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
colnames(resDE) <- paste0(res_name, c("coeffs", "t", "PValues", "AdjPValues", "Results"))

resGenes <- fit2$genes

resExpr <- round(exprs(eset_main_org), 2)
colnames(resExpr) <- paste0(targets_org$labels, "_", colnames(resExpr))
resExpr <- resExpr[, order(colnames(resExpr))]

resAll <- cbind(resGenes, resDE, resExpr)

write.table(resAll, file = paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




### plot MA
pdf(paste0(path_plots, "Comp", comp_name ,"_plotMA.pdf"))

limma::plotMA(fit2, coef = 1, status = results, values = c(-1, 0, 1), col = c("red", "black", "green"), cex = c(0.7, 0.3, 0.7))
abline(0,0,col="blue")

dev.off()




### volcano plots
library(ggplot2)


table <- topTable(fit2, coef = 1, n=Inf)
table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
gg2 <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + theme_bw() + theme(legend.position = "none") +  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("after Treatment")

pdf(paste0(path_plots, "Comp", comp_name ,"_volcanoplot.pdf"))
print(gg2)
dev.off()



### histograms of p-values and adjusted p-values
colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp", comp_name ,"_hist_pvs.pdf"))

table <- topTable(fit2, coef = 1, n=Inf)
hist(table$P.Value, breaks = 100, main = "afterTreatment", xlab = "P-values", col = colours)

dev.off()




### plot expression of ALL sign. genes/probesets
library(ggplot2)
library(reshape2)

expr <- exprs(eset_main)
rownames(targets) <- strsplit2(targets$FileName, split = "//")[, 2] 

tt <- topTable(fit2, coef = 1, n=Inf, p.value=0.05, lfc=1)

write.table(tt, paste0(path_plots, "Comp", comp_name ,"_allExpression.xls"), quote = FALSE, sep = "\t", row.names = FALSE)


topn <- nrow(tt)
top_genes <- rownames(tt)[1:topn]


df <- data.frame(Gene = top_genes, expr[top_genes,])
dfm <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")

### keep order of genes as in tt
dfm$Gene <- factor(dfm$Gene, levels = top_genes)


### add GeneSymbol and GeneTitle to the facet labels
dfm$GeneFullName <- dfm$Gene

GeneSymbol <- strsplit2(tt[top_genes, "GeneSymbol_Ensembl"], " /// ")[,1]
GeneTitle <- paste0(substr(strsplit2(tt[top_genes,"GeneTitle_Ensembl"], " /// ")[,1], 1, 30))

levels(dfm$GeneFullName) <- paste0(top_genes, " / ", GeneSymbol, "\n", GeneTitle)


dfm$groups <- factor(targets[dfm$Sample ,"groups"])

fill_colors <- unique(targets[, c("groups", "colors")])
rownames(fill_colors) <- fill_colors$groups

fill_colors <- fill_colors[levels(dfm$groups), "colors"]


Genes <- levels(dfm$Gene)
GeneFullNames <- levels(dfm$GeneFullName)

all(levels(dfm$Sample) == rownames(targets))

### Order by group
pdf(paste0(path_plots, "Comp", comp_name ,"_allExpressionBarPlot.pdf"), 4, 4)

for(i in 1:nlevels(dfm$Gene)){
  
  if(GeneSymbol[i] == "---" || grepl("immunoglobulin kappa", GeneTitle[i]) || grepl("predicted gene", GeneTitle[i]))
    next
  
  ggp <- ggplot(dfm[dfm$Gene == Genes[i], , drop = FALSE], aes(x = Sample, y = Expression, fill = groups)) +  
    geom_bar(stat = "identity") +
    theme_bw() +
    ggtitle(GeneFullNames[i]) +
    ylab("Log2 expression") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(size = 12), strip.text.x = element_text(size = 10), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank(), axis.title.y = element_text(face = "bold")) +
    scale_x_discrete(labels = targets[levels(dfm$Sample), "CellTypeShort"]) +
    scale_fill_manual(values = fill_colors)
  
  print(ggp)   
  
}

dev.off()


### Order by cell type

order_cell_type <- order(targets$CellTypeShort, targets$ExperimentShort, targets$batch)
dfm$Sample <- factor(dfm$Sample, levels = rownames(targets[order_cell_type, ]))

dfm$batch <- targets[match(dfm$Sample, rownames(targets)) ,"batch"]
dfm$batch <- factor(dfm$batch)


pdf(paste0(path_plots, "Comp", comp_name ,"_allExpressionBarPlot2.pdf"), 4, 4)

for(i in 1:nlevels(dfm$Gene)){
  
  if(GeneSymbol[i] == "---" || grepl("immunoglobulin kappa", GeneTitle[i]) || grepl("predicted gene", GeneTitle[i]))
    next
  
  ggp <- ggplot(dfm[dfm$Gene == Genes[i], , drop = FALSE], aes(x = Sample, y = Expression, fill = groups, alpha = batch)) +  
    geom_bar(stat = "identity") +
    theme_bw() +
    ggtitle(GeneFullNames[i]) +
    ylab("Log2 expression") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(size = 12), strip.text.x = element_text(size = 10), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank(), axis.title.y = element_text(face = "bold")) +
    scale_x_discrete(labels = targets[levels(dfm$Sample), "CellTypeShort"]) +
    scale_fill_manual(values = fill_colors) +
    scale_alpha_manual(values = c(1, 0.6))
  
  print(ggp)   
  
}

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

top_genes <- rownames(tt)[1:topn]

df <- data.frame(Gene = top_genes, expr[top_genes,])
dfm <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")
### keep order of genes as in tt
dfm$Gene <- factor(dfm$Gene, levels = top_genes)
### add Entrez ID to the facet labels
lab.fct <- paste0(top_genes, "\n", strsplit2(tt[top_genes, "GeneSymbol"], " /// ")[,1])
levels(dfm$Gene) <- lab.fct
dfm$groups <- targets[dfm$Sample ,"groups"]

fill_colors <- unique(targets[, c("groups", "colors")])
fill_colors <- fill_colors[order(fill_colors$groups), "colors"]

ggp <- ggplot(dfm, aes(x = Sample, y = Expression, fill = groups)) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(size = 16), strip.text.x = element_text(size = 10), legend.position = "bottom") +
  scale_x_discrete(labels=targets$CellTypeShort) +
  labs(y = "Log2 expression") +
  geom_bar(stat = "identity") +
  facet_wrap(~ Gene, scales="free_y", ncol= 5) +
  scale_fill_manual(values = fill_colors)

pdf(paste0(path_plots, "Comp", comp_name ,"_topExpressionBarPlot.pdf"), 12, 8)
print(ggp)    
dev.off()




###########################################################################
#### Comparison 3c:  pre VS after treatment with matched samples + cell type
### DO NOT use the B2M3 technical replicate from batch 2 
###########################################################################

res_name <- "afterTreatment_matched_paired_noTechnical_"
comp_name <- "3c"


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
samples2keep <- grepl("leukemia|afterTreatment", targets$labels) & targets$CellTypeShort %in% cell_types & !(targets$CellTypeShort == "B2M3" & targets$batch == "batch2" & targets$TypeOfReplicate == "technical")

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups, targets$CellTypeShort)
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


pdf(paste0(path_plots, "Comp", comp_name ,"_plotSA_trend.pdf"))
plotSA(fit2)
dev.off()


## with the FC cutoff
results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=1)
summary(results)


colours <- unique(targets[targets$groups == "afterTreatment", "colors"])

pdf(paste0(path_plots, "Comp", comp_name ,"_vennDiagram.pdf"))
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
colnames(resDE) <- paste0(res_name, c("coeffs", "t", "PValues", "AdjPValues", "Results"))

resGenes <- fit2$genes

resExpr <- round(exprs(eset_main_org), 2)
colnames(resExpr) <- paste0(targets_org$labels, "_", colnames(resExpr))
resExpr <- resExpr[, order(colnames(resExpr))]

resAll <- cbind(resGenes, resDE, resExpr)

write.table(resAll, file = paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




### plot MA
pdf(paste0(path_plots, "Comp", comp_name ,"_plotMA.pdf"))

limma::plotMA(fit2, coef = 1, status = results, values = c(-1, 0, 1), col = c("red", "black", "green"), cex = c(0.7, 0.3, 0.7))
abline(0,0,col="blue")

dev.off()




### volcano plots
library(ggplot2)


table <- topTable(fit2, coef = 1, n=Inf)
table$threshold = as.factor(table$adj.P.Val < 0.05 & abs(table$logFC) > 1)
gg2 <- ggplot(data=table, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + theme_bw() + theme(legend.position = "none") +  xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle("after Treatment")

pdf(paste0(path_plots, "Comp", comp_name ,"_volcanoplot.pdf"))
print(gg2)
dev.off()



### histograms of p-values and adjusted p-values
colours <- unique(targets[targets$groups != "leukemia", "colors"])

pdf(paste0(path_plots, "Comp", comp_name ,"_hist_pvs.pdf"))

table <- topTable(fit2, coef = 1, n=Inf)
hist(table$P.Value, breaks = 100, main = "afterTreatment", xlab = "P-values", col = colours)

dev.off()




### plot expression of ALL sign. genes/probesets
library(ggplot2)
library(reshape2)

expr <- exprs(eset_main)
rownames(targets) <- strsplit2(targets$FileName, split = "//")[, 2] 

tt <- topTable(fit2, coef = 1, n=Inf, p.value=0.05, lfc=1)

write.table(tt, paste0(path_plots, "Comp", comp_name ,"_allExpression.xls"), quote = FALSE, sep = "\t", row.names = FALSE)


topn <- nrow(tt)
top_genes <- rownames(tt)[1:topn]


df <- data.frame(Gene = top_genes, expr[top_genes,])
dfm <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")

### keep order of genes as in tt
dfm$Gene <- factor(dfm$Gene, levels = top_genes)


### add GeneSymbol and GeneTitle to the facet labels
dfm$GeneFullName <- dfm$Gene

GeneSymbol <- strsplit2(tt[top_genes, "GeneSymbol_Ensembl"], " /// ")[,1]
GeneTitle <- paste0(substr(strsplit2(tt[top_genes,"GeneTitle_Ensembl"], " /// ")[,1], 1, 30))

levels(dfm$GeneFullName) <- paste0(top_genes, " / ", GeneSymbol, "\n", GeneTitle)


dfm$groups <- factor(targets[dfm$Sample ,"groups"])

fill_colors <- unique(targets[, c("groups", "colors")])
rownames(fill_colors) <- fill_colors$groups

fill_colors <- fill_colors[levels(dfm$groups), "colors"]


Genes <- levels(dfm$Gene)
GeneFullNames <- levels(dfm$GeneFullName)

all(levels(dfm$Sample) == rownames(targets))

### Order by group
pdf(paste0(path_plots, "Comp", comp_name ,"_allExpressionBarPlot.pdf"), 4, 4)

for(i in 1:nlevels(dfm$Gene)){
  
  if(GeneSymbol[i] == "---" || grepl("immunoglobulin kappa", GeneTitle[i]) || grepl("predicted gene", GeneTitle[i]))
    next
  
  ggp <- ggplot(dfm[dfm$Gene == Genes[i], , drop = FALSE], aes(x = Sample, y = Expression, fill = groups)) +  
    geom_bar(stat = "identity") +
    theme_bw() +
    ggtitle(GeneFullNames[i]) +
    ylab("Log2 expression") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(size = 12), strip.text.x = element_text(size = 10), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank(), axis.title.y = element_text(face = "bold")) +
    scale_x_discrete(labels = targets[levels(dfm$Sample), "CellTypeShort"]) +
    scale_fill_manual(values = fill_colors)
  
  print(ggp)   
  
}

dev.off()


### Order by cell type

order_cell_type <- order(targets$CellTypeShort, targets$ExperimentShort, targets$batch)
dfm$Sample <- factor(dfm$Sample, levels = rownames(targets[order_cell_type, ]))

dfm$batch <- targets[match(dfm$Sample, rownames(targets)) ,"batch"]
dfm$batch <- factor(dfm$batch)


pdf(paste0(path_plots, "Comp", comp_name ,"_allExpressionBarPlot2.pdf"), 4, 4)

for(i in 1:nlevels(dfm$Gene)){
  
  if(GeneSymbol[i] == "---" || grepl("immunoglobulin kappa", GeneTitle[i]) || grepl("predicted gene", GeneTitle[i]))
    next
  
  ggp <- ggplot(dfm[dfm$Gene == Genes[i], , drop = FALSE], aes(x = Sample, y = Expression, fill = groups, alpha = batch)) +  
    geom_bar(stat = "identity") +
    theme_bw() +
    ggtitle(GeneFullNames[i]) +
    ylab("Log2 expression") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(size = 12), strip.text.x = element_text(size = 10), legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank(), axis.title.y = element_text(face = "bold")) +
    scale_x_discrete(labels = targets[levels(dfm$Sample), "CellTypeShort"]) +
    scale_fill_manual(values = fill_colors) +
    scale_alpha_manual(values = c(1, 0.6))
  
  print(ggp)   
  
}

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

top_genes <- rownames(tt)[1:topn]

df <- data.frame(Gene = top_genes, expr[top_genes,])
dfm <- reshape2::melt(df, id.vars = "Gene", value.name = "Expression", variable.name = "Sample")
### keep order of genes as in tt
dfm$Gene <- factor(dfm$Gene, levels = top_genes)
### add Entrez ID to the facet labels
lab.fct <- paste0(top_genes, "\n", strsplit2(tt[top_genes, "GeneSymbol"], " /// ")[,1])
levels(dfm$Gene) <- lab.fct
dfm$groups <- targets[dfm$Sample ,"groups"]

fill_colors <- unique(targets[, c("groups", "colors")])
fill_colors <- fill_colors[order(fill_colors$groups), "colors"]

ggp <- ggplot(dfm, aes(x = Sample, y = Expression, fill = groups)) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), plot.title = element_text(size = 16), strip.text.x = element_text(size = 10), legend.position = "bottom") +
  scale_x_discrete(labels=targets$CellTypeShort) +
  labs(y = "Log2 expression") +
  geom_bar(stat = "identity") +
  facet_wrap(~ Gene, scales="free_y", ncol= 5) +
  scale_fill_manual(values = fill_colors)

pdf(paste0(path_plots, "Comp", comp_name ,"_topExpressionBarPlot.pdf"), 12, 8)
print(ggp)    
dev.off()






###########################################################################
#### Merge all results
###########################################################################


res_files <- c("Comp3b_DE_results_All.xls", "Comp3c_DE_results_All.xls")


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
vennDiagram(results, include = c("up", "down"), counts.col = c("gold", "darkblue"), cex=c(0.7,1,0.7))
dev.off()










###########################################################################
#### GO analysis for comparison 3c
###########################################################################

res_name <- "afterTreatment_matched_paired_noTechnical_"
comp_name <- "3c"


allLines <- readLines(paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), n = -1)[-1]
resAll <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(resAll) <- strsplit2(readLines(paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), n = 1), "\t")



affyLib <- "mogene20sttranscriptcluster.db"

ls("package:mogene20sttranscriptcluster.db")


### Function used to create new topGOdata object
fun.gene.sel <- function(geneList) {
  geneList <- ifelse(geneList == 0, FALSE, TRUE)
  return(geneList)
}


tt <- resAll[, c("ProbesetID", paste0(res_name, "Results"))]
geneList <- as.numeric(as.numeric(tt[, 2]) != 0)
names(geneList) <- tt[, 1]
print(table(geneList))

table(fun.gene.sel(geneList))

geneList_sign <- names(geneList[geneList != 0])


for(go in c("BP","MF","CC")){
  # go = "BP"
  print(go)
  
  sampleGOdata <- new("topGOdata", description = paste0("Simple session for ", comp_name), ontology = go, allGenes = geneList, geneSel = fun.gene.sel , nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
  
  print(sampleGOdata)
  
  result <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
  
  pValues <- score(result)
  topNodes <- length(pValues)
  
  allRes <- GenTable(sampleGOdata, elimFisher = result, orderBy = "elimFisher", topNodes = topNodes, numChar = 300) 
  
  colnames(allRes)[6] <- "PValues" 
  
  allRes$PValues <- as.numeric(allRes$PValues)
  allRes$Significant <- as.numeric(allRes$Significant)
  
  
  #   pdf(paste(path_plots, "GO_", go, ".pdf", sep=""))
  #   showSigOfNodes(sampleGOdata, score(result), firstSigNodes = 5, useInfo = 'all')
  #   dev.off()
  #   
  #   allRes$AdjPValues <- p.adjust(allRes$PValues, method = "BH")
  #   
  #   print(table(allRes$AdjPValues < 0.05))
  #   
  #   cat("#########################################################################################", fill = TRUE)
  #   print(head(allRes, 20))
  #   cat("#########################################################################################", fill = TRUE)
  
  
  which_go <- which(allRes$PValues < 0.05 & allRes$Significant >= 10)
  
  genesInNode <- sapply(which_go, function(i){
    # i = 1
    goID <- allRes[i, "GO.ID"]
    goID
    
    # gt <- printGenes(sampleGOdata, whichTerms = goID, chip = affyLib)[, 1]
    
    gt <- genesInTerm(sampleGOdata, goID)[[1]]
    
    paste0(as.character(gt[gt %in% geneList_sign]), collapse = " /// ")
    
    
  })
  
  
  allRes$genesInNode <- "---"
  allRes$genesInNode[which_go] <- genesInNode
  
  
  #   termStat(sampleGOdata, whichGO = goID)
  #   
  #   pdf(paste(path_plots, "GO_", go, ".pdf", sep=""))
  #   showGroupDensity(sampleGOdata, whichGO = goID)
  #   dev.off()
  
  
  write.table(allRes, paste0(path_results, "Comp", comp_name ,"_GO_Fisher_elim_", go ,".xls"), sep="\t", row.names=FALSE, quote = FALSE)
  
  
}


###########################################################################
#### Add GO terms to the significant genes in comparison 3c based genesInTerm - GO graph
###########################################################################


res_name <- "afterTreatment_matched_paired_noTechnical_"
comp_name <- "3c"


allLines <- readLines(paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), n = -1)[-1]
resAll <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(resAll) <- strsplit2(readLines(paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), n = 1), "\t")


resSign <- resAll[as.numeric(resAll[, paste0(res_name, "Results")]) != 0,  !grepl(".CEL", colnames(resAll))]



for(go in c("BP","MF","CC")){
  # go = "CC"
  
  allRes <- read.table(paste0(path_results, "Comp", comp_name ,"_GO_Fisher_elim_", go ,".xls"), header = TRUE, sep = "\t", as.is = TRUE)
  
  allRes <- allRes[allRes$genesInNode != "---", ]
  
  genesInNode_split <- strsplit(allRes$genesInNode, split = " /// ")
  
  names(genesInNode_split) <- allRes$Term
  
  
  mysets_df <- data.frame(GO_terms = rep(names(genesInNode_split), sapply(genesInNode_split, length)), ProbesetID = unlist(genesInNode_split), row.names = NULL, stringsAsFactors = FALSE)
  
  
  GO_terms <- sapply(1:nrow(resSign), function(i){
    # i = 1
    
    if(!resSign$ProbesetID[i] %in% mysets_df$ProbesetID)
      return("---")
    
    
    TermsDefinition <- mysets_df[mysets_df$ProbesetID == resSign$ProbesetID[i], "GO_terms"]
    
    paste0(unique(TermsDefinition), collapse = " /// ")
    
    
  })
  
  
  
  resSign[, paste0("GO_terms_", go)] <- GO_terms
  
  
}


write.table(resSign, file = paste0(path_results, "Comp", comp_name ,"_DE_results_Sign.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




###########################################################################
#### Add GO terms to the significant genes in comparison 3c based on annotation 
### gives less genes than shown in allRes
###########################################################################

res_name <- "afterTreatment_matched_paired_noTechnical_"
comp_name <- "3c"


allLines <- readLines(paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), n = -1)[-1]
resAll <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(resAll) <- strsplit2(readLines(paste0(path_results, "Comp", comp_name ,"_DE_results_All.xls"), n = 1), "\t")


resSign <- resAll[as.numeric(resAll[, paste0(res_name, "Results")]) != 0,  !grepl(".CEL", colnames(resAll))]


### get GO terms from MSigDB
# gene sets from MSigDB with ENTREZ IDs
load("MSigDB_v4_0/mouse_c5_v4.rdata")

mysets <- Mm.c5
length(mysets)


mysets_df <- data.frame(MSigDB_GO_terms = rep(names(mysets), sapply(mysets, length)), EntrezGeneID = unlist(mysets), row.names = NULL, stringsAsFactors = FALSE)

EntrezGeneID <- strsplit2(resSign$EntrezGeneID, " /// ")


MSigDB_GO_terms <- sapply(1:nrow(resSign), function(i){
  # i = 2
  
  if(EntrezGeneID[i, 1] == "---")
    return("---")
  
  TermsDefinition <- mysets_df$MSigDB_GO_terms[which(mysets_df$EntrezGeneID %in% EntrezGeneID[i, ])]
  
  if(length(TermsDefinition) == 0)
    return("---")
  
  paste0(unique(TermsDefinition), collapse = " /// ")
  
  
})


# resSign$MSigDB_GO_terms <- MSigDB_GO_terms
# 
# 
# write.table(resSign, file = paste0(path_results, "Comp", comp_name ,"_DE_results_Sign.xls"), quote = FALSE, sep = "\t", row.names = FALSE)




### get GO terms from mogene20sttranscriptcluster.db


ls("package:mogene20sttranscriptcluster.db")

x <- mogene20sttranscriptclusterGO
# Get the manufacturer identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)

all(resSign$ProbesetID %in% mapped_genes)

# Convert to a list
xx <- as.list(x[mapped_genes])


GO_terms <- sapply(1:nrow(resSign), function(i){
  # i = 28
  
  if(!resSign$ProbesetID[i] %in% mapped_genes)
    return("---")
  
  whichTerms <- names(xx[[resSign$ProbesetID[i]]])
  
  TermsDefinition <- c(topGO:::.getTermsDefinition(whichTerms, ontology = "BP", numChar = 300, multipLines = FALSE),
                       topGO:::.getTermsDefinition(whichTerms, ontology = "CC", numChar = 300, multipLines = FALSE),
                       topGO:::.getTermsDefinition(whichTerms, ontology = "MF", numChar = 300, multipLines = FALSE))
  
  TermsDefinition <- TermsDefinition[!is.na(TermsDefinition)]
  
  paste0(unique(TermsDefinition), collapse = " /// ")
  
  
})


# resSign$GO_terms <- GO_terms
# 
# 
# write.table(resSign, file = paste0(path_results, "Comp", comp_name ,"_DE_results_Sign.xls"), quote = FALSE, sep = "\t", row.names = FALSE)





### get GO terms from mogene20sttranscriptcluster.db + child terms


x <- mogene20sttranscriptclusterGO2ALLPROBES

mapped_gos <- mappedkeys(x)

# Convert to a list
xx <- as.list(x[mapped_gos])


mysets_df <- data.frame(GO_terms = rep(names(xx), sapply(xx, length)), ProbesetID = unlist(xx), stringsAsFactors = FALSE)


mysets_df_tmp <- mysets_df[mysets_df$GO_terms == "GO:0002455", ]





GO_terms <- sapply(1:nrow(resSign), function(i){
  # i = 28
  
  if(!resSign$ProbesetID[i] %in% mysets_df$ProbesetID)
    return("---")
  
  
  whichTerms <- mysets_df$GO_terms[which(mysets_df$ProbesetID == resSign$ProbesetID[i])]
  
  TermsDefinition <- c(topGO:::.getTermsDefinition(whichTerms, ontology = "BP", numChar = 300, multipLines = TRUE),
                       topGO:::.getTermsDefinition(whichTerms, ontology = "CC", numChar = 300, multipLines = TRUE),
                       topGO:::.getTermsDefinition(whichTerms, ontology = "MF", numChar = 300, multipLines = TRUE))
  
  TermsDefinition <- TermsDefinition[!is.na(TermsDefinition)]
  
  paste0(unique(TermsDefinition), collapse = " /// ")
  
  
})



GO_ids <- sapply(1:nrow(resSign), function(i){
  # i = 28
  
  if(!resSign$ProbesetID[i] %in% mysets_df$ProbesetID)
    return("---")
  
  whichTerms <- mysets_df$GO_terms[which(mysets_df$ProbesetID == resSign$ProbesetID[i])]
  
  paste0(unique(whichTerms), collapse = " /// ")
  
  
})




# 
# resSign$GO_terms_child <- GO_terms
# resSign$GO_ids_child <- GO_ids
# 
# 
# write.table(resSign, file = paste0(path_results, "Comp", comp_name ,"_DE_results_Sign.xls"), quote = FALSE, sep = "\t", row.names = FALSE)












###########################################################################
#### Gene set enrichment analysis with CAMERA for C5 - GO genes sets
# http://bioinf.wehi.edu.au/software/MSigDB/
### comparison 3c
###########################################################################

res_name <- "afterTreatment_matched_paired_noTechnical_"
comp_name <- "3c"

load(paste0(path_results, "eset_main_org.Rdata"))
load(paste0(path_results, "targets_org.Rdata"))

targets <- targets_org
eset_main <- eset_main_org

tt <- table(targets$CellTypeShort, targets$groups)

tt <- data.frame(tt, stringsAsFactors = FALSE)

cell_types <- as.character(tt[tt$Var2 == "afterTreatment" & tt$Freq > 0, "Var1"])

### keep only leukemia and afterTreatment samples that have matched cell type 
samples2keep <- grepl("leukemia|afterTreatment", targets$labels) & targets$CellTypeShort %in% cell_types & !(targets$CellTypeShort == "B2M3" & targets$batch == "batch2" & targets$TypeOfReplicate == "technical")

targets <- targets[samples2keep,]
eset_main <- eset_main[, samples2keep]

### sort samples by groups
ord <- order(targets$groups, targets$CellTypeShort)
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


EntrezGeneID <- strsplit2(annot$EntrezGeneID, " /// ")

nrow <- nrow(EntrezGeneID)
ncol <- ncol(EntrezGeneID)

Index <- lapply(mysets, function(ms){  
  eg <- matrix(EntrezGeneID %in% ms, nrow = nrow, ncol = ncol, byrow = FALSE)
  rowSums(eg) > 0 
})




### run CAMERA

gsea <- camera(y = eset_main, index = Index, design = design, contrast = ncol(design), trend.var = TRUE)
head(gsea, 10)
table(gsea$FDR < 0.05)

gsea <- data.frame(GeneSet = rownames(gsea), gsea)

write.table(gsea, paste(path_results, "Comp", comp_name ,"_GSEA_C5_GO_gene_sets.xls", sep=""), sep="\t", row.names=F, quote = FALSE)





















