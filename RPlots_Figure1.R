###########################################################################
# Created 20 Apr 2015
# BioC 3.0

# Generate plots for Figure 1

###########################################################################

setwd("/home/Shared/data/array/Microarray_Edwin")

###########################################################################
### data 
###########################################################################

targets <- read.table(file.path("metadata", "targets.xls"), header = TRUE, sep = "\t", comment.char = "", as.is = TRUE)

targets

library(oligo)
library(pd.mogene.2.0.st)

ff <- as.character(targets$FileName)

x <- oligo::read.celfiles(filenames = ff)


# Normalize data with RMA method

eset <- oligo::rma(x) 

####### Keep only the control and leukemia samples

keepSAMPS <- targets$ExperimentShort != "afterTreatment" & targets$labels != "control_HeLa"

eset <- eset[, keepSAMPS]
targets <- targets[keepSAMPS, ]


###########################################################################
### MDS plot
###########################################################################

library(limma)
labels <- paste0(ifelse(is.na(targets$ctrlRep), "", paste0(targets$ctrlRep, " ")), targets$labels)

### plot MDS for all samples
mds <- plotMDS(eset, top=1000, col = targets$colors, labels = labels, cex = 1.2)


###########################################################################
### Heatmaps
###########################################################################

library(pheatmap)


### sort samples by groups
ord <- order(targets$groups)
targets <- targets[ord, ]
eset <- eset[ ,ord]



resAll <- read.table("Comp1_DE_results_All.xls", header = TRUE, sep = "\t", stringsAsFactors = FALSE, nrows = 5)


allLines <- readLines("Comp1_DE_results_All.xls", n = -1)[-1]
resAll <- data.frame(strsplit2(allLines, "\t"), stringsAsFactors = FALSE)
colnames(resAll) <- strsplit2(readLines("Comp1_DE_results_All.xls", n = 1), "\t")

expr <- data.frame(exprs(eset))

### load list of genes of interest

intrGenes <- read.table("Gene_Sets/niche_interaction_genes_niche_binding.csv", header = TRUE, sep = ",")
intrProbes <- as.character(intrGenes$ProbesetID)

dataHeat <- expr[intrProbes, ]

annotation_col <- targets[, "groups", drop = FALSE]
rownames(annotation_col) <- colnames(dataHeat)

cols <- unique(targets$colors)
names(cols) <- unique(targets$group)
  
annotation_colors = list(groups = cols)

labels_row <- strsplit2(intrGenes$GeneSymbol, " /// ")[, 1]



pdf("PLOTS/heatmap.pdf", width = 7, height = 10)

pheatmap(dataHeat, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 7, breaks =  seq(from = 2, to = 14, length.out = 101), legend_breaks = seq(from = 2, to = 14, by = 2))

pheatmap(dataHeat, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 7, scale = "row")

dev.off()



### load list of genes of interest

intrGenes <- read.table("Gene_Sets/niche_interaction_genes_niche_formation.csv", header = TRUE, sep = ",")
intrProbes <- as.character(intrGenes$ProbesetID)

dataHeat <- expr[intrProbes, ]

annotation_col <- targets[, "groups", drop = FALSE]
rownames(annotation_col) <- colnames(dataHeat)

cols <- unique(targets$colors)
names(cols) <- unique(targets$group)

annotation_colors = list(groups = cols)

labels_row <- strsplit2(intrGenes$GeneSymbol, " /// ")[, 1]



pdf("PLOTS/heatmap.pdf", width = 7, height = 10)

pheatmap(dataHeat, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 7, breaks =  seq(from = 2, to = 14, length.out = 101), legend_breaks = seq(from = 2, to = 14, by = 2))

pheatmap(dataHeat, cluster_cols = FALSE, annotation_col = annotation_col, annotation_colors = annotation_colors, labels_col = targets$groups, labels_row = labels_row, annotation_legend = FALSE, fontsize_row = 7, scale = "row")


dev.off()
















































