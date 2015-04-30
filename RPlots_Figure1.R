###########################################################################
# Created 20 Apr 2015
# BioC 3.0

# Generate plots for Figure 1

###########################################################################

setwd("/home/Shared/data/array/Microarray_Edwin")

###########################################################################
### data 
###########################################################################

targets.org <- read.table(file.path("metadata", "targets.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)

targets.org

library(oligo)
library(pd.mogene.2.0.st)

ff <- as.character(targets.org$FileName)

x <- oligo::read.celfiles(filenames = ff)


# Normalize data with RMA method

eset.org <- oligo::rma(x) 

####### Keep only the control and leukemia samples

keepSAMPS <- targets.org$ExperimentShort != "afterTreatment" & targets.org$labels != "control_HeLa" & targets.org$labels != "control_wholeBoneMarrow"

eset <- eset.org[, keepSAMPS]
targets <- targets.org[keepSAMPS, ]


###########################################################################
### MDS plot
###########################################################################

library(limma)
library(ggplot2)


labels <- paste0(ifelse(is.na(targets$ctrlRep), "", paste0(targets$ctrlRep, " ")), targets$labels)



### plot MDS for all samples
pdf("PLOTS_Fig1/MDS.pdf", width = 7, height = 7)
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.5, xlim = c(-2, 3), ylim = c(-2, 3), las = 1, cex.axis = 1.5, cex.lab = 1.5)
dev.off()


pdf("PLOTS_Fig1/MDS.pdf", width = 5, height = 5)

mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1, xlim = c(-2, 3), ylim = c(-2, 3), las = 1, cex.axis = 1, cex.lab = 1)

dev.off()


# df <- data.frame(x = mds$x, y = mds$y, labels = factor(labels), colors = factor(targets$colors))
# 
# 
# ggp2 <- ggplot(df, aes(x = x, y = y, colour = colors)) +
#   theme_bw() +
#   geom_point(size = 2) +
#   geom_text(aes(label=labels), hjust = 0.5, vjust = 0.5) +
#   theme(legend.position="none")+
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 12), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
#   labs(x = "Leading logFC dim 1", y = "Leading logFC dim 2") +
#   coord_cartesian(xlim = c(-2, 3), ylim = c(-2, 3)) +
#   scale_colour_manual(values = levels(df$colors))
# 
# 
# 
# pdf("PLOTS_Fig1/MDSgg.pdf", width = 5, height = 5)
# 
# print(ggp2)
# 
# dev.off()



######### MDS plot based on DE genes

allLines <- readLines("Comp1_DE_results_All.xls", n = -1)[-1]
resAll <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(resAll) <- strsplit2(readLines("Comp1_DE_results_All.xls", n = 1), "\t")

resAll <- resAll[, !grepl(pattern = "CEL", colnames(resAll))]

resAllSort <- resAll[order(resAll$CtrlCD4_PValues, resAll$CtrlCD4CD8_PValues, resAll$CtrlCD8_PValues, decreasing = FALSE), ]

resAllSort$clusters <- apply(resAllSort[, c("CtrlCD4_Results", "CtrlCD4CD8_Results", "CtrlCD8_Results")], MARGIN = 1, paste, collapse = ",")

resAllSort <- resAllSort[resAllSort$clusters != "0,0,0", ]




pdf("PLOTS_Fig1/MDS_DE.pdf", width = 5, height = 5)

mds <- plotMDS(eset[resAllSort$ProbesetID, ], top=nrow(resAllSort), col = targets$colors, labels = labels, cex = 1, xlim = c(-2, 2), ylim = c(-2, 2), las = 1, cex.axis = 1, cex.lab = 1)

dev.off()





###########################################################################
### Heatmaps
###########################################################################

library(pheatmap)
library(limma)

### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset <- eset[ ,ord]


expr <- data.frame(exprs(eset))
### normalize expression per gene
exprNorm <- t(scale(t(expr), center = TRUE, scale = TRUE))

### load list of genes of interest

intrGenes <- read.table("Gene_Sets/niche_interaction_genes_v3_niche_binding.csv", header = TRUE, sep = ",")
intrGenes <- read.table("Gene_Sets/niche_interaction_genes_v3_niche_formation.csv", header = TRUE, sep = ",")
intrProbes <- as.character(intrGenes$ProbesetID)

# dataHeat <- expr[intrProbes, ]
dataHeat <- exprNorm[intrProbes, ]



annotation_col <- targets[, "groups", drop = FALSE]
colnames(annotation_col) <- "Samps"
rownames(annotation_col) <- colnames(dataHeat)

cols <- unique(targets$colors)
names(cols) <- unique(targets$group)
  
annotation_colors = list(groups = cols)
names(annotation_colors) <- "Samps"


labels_row <- strsplit2(intrGenes$GeneSymbol, " /// ")[, 1]



pdf("PLOTS_Fig1/heatmap.pdf", width = 5.5, height = 5.5)

pheatmap(dataHeat, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 6, fontsize_col = 6, gaps_col = cumsum(table(targets$groups)), breaks =  seq(from = -3, to = 3, length.out = 101), legend_breaks = seq(from = -3, to = 3, by = 2), main = "Niche binding", cellheight = 7, cellwidth = 14, fontsize = 7)

dev.off()



library(RColorBrewer)
library(colorRamps)

pdf("PLOTS_Fig1/heatmap.pdf", width = 5, height = 5)

pheatmap(dataHeat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 6, fontsize_col = 6, gaps_col = cumsum(table(targets$groups)), breaks =  seq(from = -3, to = 3, length.out = 101), legend_breaks = seq(from = -3, to = 3, by = 2), main = "Niche binding", fontsize = 7)

pheatmap(dataHeat, color = colorRamps::matlab.like(100), cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 6, fontsize_col = 6, gaps_col = cumsum(table(targets$groups)), breaks =  seq(from = -3, to = 3, length.out = 101), legend_breaks = seq(from = -3, to = 3, by = 2), main = "Niche binding", fontsize = 7)


pheatmap(dataHeat, color = colorRamps::matlab.like2(100), cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 6, fontsize_col = 6, gaps_col = cumsum(table(targets$groups)), breaks =  seq(from = -3, to = 3, length.out = 101), legend_breaks = seq(from = -3, to = 3, by = 2), main = "Niche binding", fontsize = 7)



dev.off()





































