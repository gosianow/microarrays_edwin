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
labels <- paste0(ifelse(is.na(targets$ctrlRep), "", paste0(targets$ctrlRep, " ")), targets$labels)

### plot MDS for all samples
pdf("PLOTS_Fig1/MDS.pdf")
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2, xlim = c(-2, 3))
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
exprNorm2 <- t(scale(t(expr), center = TRUE, scale = FALSE))


### load list of genes of interest

intrGenes <- read.table("Gene_Sets/niche_interaction_genes_v2_niche_binding.csv", header = TRUE, sep = ",")
intrProbes <- as.character(intrGenes$ProbesetID)

# dataHeat <- expr[intrProbes, ]
dataHeat <- exprNorm[intrProbes, ]

annotation_col <- targets[, "groups", drop = FALSE]
rownames(annotation_col) <- colnames(dataHeat)

cols <- unique(targets$colors)
names(cols) <- unique(targets$group)
  
annotation_colors = list(groups = cols)

labels_row <- strsplit2(intrGenes$GeneSymbol, " /// ")[, 1]



pdf("PLOTS_Fig1/heatmap.pdf", width = 7, height = 10)

pheatmap(dataHeat, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 8, gaps_col = cumsum(table(targets$groups)), breaks =  seq(from = -3, to = 3, length.out = 101), legend_breaks = seq(from = -3, to = 3, by = 2), main = "Niche binding")

dev.off()





###########################################################################
### Clustering expression data for all genes 
###########################################################################

## does not work 
# resAll <- read.table("Comp1_DE_results_All.xls", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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

unique(resAllSort$clusters)





intrProbes <- as.character(resAllSort$ProbesetID)

# dataHeat <- expr[intrProbes, ]
dataHeat <- exprNorm[intrProbes, ]

annotation_col <- targets[, "groups", drop = FALSE]
rownames(annotation_col) <- colnames(dataHeat)

cols <- unique(targets$colors)
names(cols) <- unique(targets$group)

annotation_colors = list(groups = cols)

labels_row <- strsplit2(intrGenes$GeneSymbol, " /// ")[, 1]



png("PLOTS_Fig1/heatmap_all.png", width = 700, height = 2000)

pheatmap(dataHeat, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = rep("", nrow(dataHeat)), annotation_legend = FALSE, fontsize_row = 8, gaps_col = cumsum(table(targets$groups)), gaps_row = cumsum(table(resAllSort$clusters)),breaks =  seq(from = -4, to = 4, length.out = 101), legend_breaks = seq(from = -4, to = 4, by = 2))

dev.off()



write.table(resAllSort, file = "Comp1_DEclusters.xls", quote = FALSE, sep = "\t", row.names = FALSE)



































