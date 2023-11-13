# ================================================================
# Publication title: Polymer-mediated oligonucleotide delivery
# enables construction of spatially encoded 3D cultures for
# single-cell RNA sequencing analysis 
#
# Code authors: Jessica J. King and Cameron W. Evans
# Date: October 2023
#
# Contents:
# line 17   Preliminary declarations and setup
# line 46   Figure 1: Untreated 2-layer spheroids
# line 267  Figure 2: Untreated 3-layer spheroids
# line 543  Figure 4/5: Treated 3-layer spheroids
# ================================================================


### Load packages

library(AnnotationHub)
library(AnnotationDbi)
library(circlize)
library(clusterProfiler)
library(ComplexHeatmap)
library(dittoSeq)
library(dplyr)
library(escape)
library(fgsea)
library(ggplot2)
library(GSEABase)
library(msigdbr)
library(org.Hs.eg.db)
library(patchwork)
library(ReactomePA)
library(RColorBrewer)
library(Seurat)
library(SingleCellExperiment)


### Initial setup

currentDir <- "/Users/00058955/Library/CloudStorage/OneDrive-SharedLibraries-TheUniversityofWesternAustralia/RES-BioNano - General/Papers/Jess King - spheroids/R code/"
setwd(currentDir)



### ==============================================================
### Figure 1: Untreated 2-layer spheroids
### ==============================================================


# ----------------------------------------------------------------
# Load the data
# ----------------------------------------------------------------

load("SCdata_2layer.RData")


# ----------------------------------------------------------------
# Subset and cluster the cells
# ----------------------------------------------------------------

# Subset by sample
SCdata <- SetIdent(SCdata, value = SCdata$Sample)
SCdata <- subset(x = SCdata, idents = c("unsorted-core","unsorted-periphery"))

# Clustering
normdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
normdata <- FindVariableFeatures(normdata, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(normdata)
scaleddata <- ScaleData(normdata)
dimred <- RunPCA(scaleddata, features = VariableFeatures(object = scaleddata))
SCdata <- FindNeighbors(dimred, dims = 1:20)
SCdata <- FindClusters(SCdata, resolution = 0.1)
SCdata <- RunUMAP(SCdata, dims = 1:20)

# Clean up
rm(normdata, all.genes, scaleddata, dimred)


# ----------------------------------------------------------------
# Analyse cell cycle
# ----------------------------------------------------------------

# Scoring cell cycle
s.genes <- cc.genes.updated.2019$s.genes 
s.genes <- s.genes[!is.na(match(s.genes, rownames(SCdata)))]
g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- g2m.genes[!is.na(match(g2m.genes, rownames(SCdata)))] 
cellcycle = CellCycleScoring(SCdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 
SCdata$cellcycle <- paste(cellcycle$Phase)

# Clean up
rm(s.genes, g2m.genes, cellcycle)


# ----------------------------------------------------------------
# Check metadata
# ----------------------------------------------------------------

# Subset by sample
SCdata <- SetIdent(SCdata, value = SCdata$Sample)
SCdata_core <- subset(x = SCdata, idents = "unsorted-core")
SCdata_periphery <- subset(x = SCdata, idents = "unsorted-periphery")

# Metadata summary - transcript counts
summary(SCdata_core@meta.data[["nCount_RNA"]]) 
summary(SCdata_periphery@meta.data[["nCount_RNA"]]) 
summary(SCdata@meta.data[["nCount_RNA"]])

# Metadata summary - gene counts
summary(SCdata_core@meta.data[["nFeature_RNA"]]) 
summary(SCdata_periphery@meta.data[["nFeature_RNA"]]) 
summary(SCdata@meta.data[["nFeature_RNA"]]) 

# Metadata summary - SBO counts
summary(SCdata_core@meta.data[["nCount_HTO"]]) 
summary(SCdata_periphery@meta.data[["nCount_HTO"]]) 
summary(SCdata@meta.data[["nCount_HTO"]]) 

# Clean up
rm(SCdata_core, SCdata_periphery)


# ----------------------------------------------------------------
# Plot SBO reads (Cy5 vs Cy3) to show effective demultiplexing
# ----------------------------------------------------------------

outputFilename <- "Fig1g.csv"
f <- FeatureScatter(SCdata, feature1 = "SBO01-cy5", feature2 = "SBO02-cy3") +
  xlim(-100, 30000) + ylim(-100, 30000) + ylab("Cy3-SBO") + xlab("Cy5-SBO") +
  geom_abline(intercept = 0, slope = 1, linetype = 2) + theme(text = element_text(size = 12))
show(f)
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Plot UMI and SBO counts for Cy3 vs Cy5 to check no difference
# ----------------------------------------------------------------

# Plot distribution of UMI counts
outputFilename <- "Fig1h.csv"
f <- VlnPlot(SCdata, features = "nCount_RNA",
             pt.size = 0.1, log = TRUE, group.by = "Sample", ncol = 1) + ylab("Number of UMIs per cell") +
  theme(plot.title = element_blank(), legend.position = "none",  text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, vjust = -1))
show(f)
write.csv(f[[1]]$data, paste(currentDir,outputFilename,sep = ""))

# Plot distribution of SBO counts
outputFilename <- "Fig1i.csv"
f <- VlnPlot(SCdata, features = "nCount_HTO",
             pt.size = 0.1, log = TRUE, group.by = "Sample", ncol = 1) + ylab("SBO reads per cell") +
  theme(plot.title = element_blank(), legend.position = "none",  text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, vjust = -1))
show(f)
write.csv(f[[1]]$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# View UMAP clustering and how that compares across samples
# ----------------------------------------------------------------

# Make UMAP by cluster
f <- DimPlot(SCdata, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = F) + 
  theme(text = element_text(size = 12)) 
show(f)
outputFilename <- "2layer_UMAP_cluster.csv"
write.csv(f[[1]]$data, paste(currentDir,outputFilename,sep = ""))

# Make UMAP by sample
f <- DimPlot(SCdata, reduction = "umap", group.by = "Sample", pt.size = 0.5, label = F) + 
  theme(text = element_text(size = 12)) 
show(f)
outputFilename <- "2layer_UMAP_sample.csv"
write.csv(f[[1]]$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Make bar plots for cell cycle distribution
# ----------------------------------------------------------------

# Stacked bar plot showing cell cycle vs cluster
f <- dittoBarPlot(SCdata, "cellcycle", group.by = "seurat_clusters")
show(f)

# Stacked bar plot showing cell cycle vs layer
f <- dittoBarPlot(SCdata, "cellcycle", group.by = "Sample")
show(f)


# ----------------------------------------------------------------
# Identify markers of cluster 2
# ----------------------------------------------------------------

# Have a look and see which genes characterise cluster 2
SCdata <- SetIdent(SCdata, value = "seurat_clusters")
cluster2markers <- FindMarkers(SCdata, ident.1 = 2, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
cluster2markers <- cbind("gene" = rownames(cluster2markers),cluster2markers)
allDEGs <- unique(cluster2markers$gene)
print(paste("There are",length(allDEGs),"DEGs that identify cluster 2."))
# Conclusion: lots of genes, just about everything because the cells are apoptotic?
# Note cluster 2 is spread across both interior and periphery so doesn't impact DEGs between layers

# Clean up
rm(cluster2markers, allDEGs)


# ----------------------------------------------------------------
# Identify DEGs across samples based on average expression
# ----------------------------------------------------------------

# Remove cluster 2 if desired (doesn't influence this result)
# SCdata <- SetIdent(SCdata, value = "seurat_clusters")
# SCdata <- subset(SCdata, ident = 2, invert = TRUE)

# First get all DEGs (note doesn't seem to matter if we include or exclude cluster 2 here)
SCdata <- SetIdent(SCdata, value = "Sample")

# Find all DEGs and filter out any that have adjusted p-values > 0.05
all.markers_markers <- FindAllMarkers(SCdata, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)

allDEGs <- unique(all.markers_markers$gene)
print(paste("There are",length(allDEGs),"DEGs between the different layers."))


# ----------------------------------------------------------------
# Identify upregulated pathways in core and periphery
# ----------------------------------------------------------------

all.markers_markers <- filter(all.markers_markers, avg_log2FC > 0)
dfsample <- split(all.markers_markers$gene,all.markers_markers$cluster)

# Convert genes to ENTREZIDs
genelist <- c()
group_names <- names(dfsample)
for (x in group_names) {
  dfsample[[x]] = bitr(dfsample[[x]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genelist[[x]] <- dfsample[[x]]$ENTREZID
}
names(genelist) <- c("Interior", "Periphery")
rm(x,dfsample,group_names)

# Examine the most affected pathways
test_comparison <- compareCluster(genelist, fun = "enrichPathway")
f <- dotplot(test_comparison, x = "Cluster", color = "p.adjust", showCategory = 20,
             by = "geneRatio",
             size = "Count",
             split = NULL,
             includeAll = TRUE,
             font.size = 10,
             title = "",
             label_format = 50,
             group = FALSE,
             shape = FALSE)
show(f)

# Clean up
rm(all.markers_markers, allDEGs, genelist, test_comparison, dfsample)


# ----------------------------------------------------------------
# End of Figure 1, can delete things now
# ----------------------------------------------------------------

# Clear subsetted data
rm(SCdata_core, SCdata_periphery)

# Clear figure and output variables
rm(f, outputFilename)

# Clear SC experiment
rm(SCdata, Avg_exp)



# ================================================================
# Figure 2: Untreated 3-layer spheroids
# ================================================================


# ----------------------------------------------------------------
# Load the data
# ----------------------------------------------------------------

load("SCdata_drug_and_control_sbos.RData")


# ----------------------------------------------------------------
# Subset the data (only looking at control) for this entire figure
# ----------------------------------------------------------------

SCdata <- SetIdent(SCdata, value = SCdata$Sample)
SCdata <- subset(x = SCdata, idents = c("control-core","control-middle","control-periphery"))


# ----------------------------------------------------------------
# Process the data
# ----------------------------------------------------------------

# Clustering
normdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
normdata <- FindVariableFeatures(normdata, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(normdata)
scaleddata <- ScaleData(normdata)
dimred <- RunPCA(scaleddata, features = VariableFeatures(object = scaleddata))
SCdata <- FindNeighbors(dimred, dims = 1:20)
SCdata <- FindClusters(SCdata, resolution = 0.2)
SCdata <- RunUMAP(SCdata, dims = 1:20)

# Generate average expression by layer
Avg_exp <- AverageExpression(SCdata, assays = "RNA", return.seurat = TRUE, group.by = "Sample", slot = "data")

# Clean up
rm(normdata, all.genes, scaleddata, dimred)

# Scoring cell cycle
s.genes <- cc.genes.updated.2019$s.genes 
s.genes <- s.genes[!is.na(match(s.genes, rownames(SCdata)))]
g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- g2m.genes[!is.na(match(g2m.genes, rownames(SCdata)))] 
cellcycle = CellCycleScoring(SCdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 
SCdata$cellcycle <- paste(cellcycle$Phase)

# Clean up
rm(s.genes, g2m.genes, cellcycle)


# ----------------------------------------------------------------
# Check metadata
# ----------------------------------------------------------------

# Subset by sample
SCdata_core <- subset(x = SCdata, idents = "control-core")
SCdata_middle <- subset(x = SCdata, idents = "control-middle")
SCdata_periphery <- subset(x = SCdata, idents = "control-periphery")

# Metadata summary
summary(SCdata_core@meta.data[["nCount_RNA"]]) 
summary(SCdata_middle@meta.data[["nCount_RNA"]]) 
summary(SCdata_periphery@meta.data[["nCount_RNA"]]) 
summary(SCdata_core@meta.data[["nCount_HTO"]]) 
summary(SCdata_middle@meta.data[["nCount_HTO"]]) 
summary(SCdata_periphery@meta.data[["nCount_HTO"]]) 

# Clean up
rm(SCdata_core, SCdata_middle, SCdata_periphery)

# Check cell number in cluster 3
SCdata <- SetIdent(SCdata, value = SCdata$seurat_clusters)
cluster3pc <- length(subset(SCdata, ident = 3)$seurat_clusters)/length(subset(SCdata, ident = c(0,1,2,3))$seurat_clusters)
print(paste("Cluster 3 represents",round(cluster3pc*100,2),"% of all cells."))


# ----------------------------------------------------------------
# Prepare UMAP plots
# ----------------------------------------------------------------

# Make UMAP by cluster
f <- DimPlot(SCdata, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = F) + 
  theme(text = element_text(size = 12)) 
show(f)
outputFilename <- "Control_UMAP_cluster.csv"
write.csv(f[[1]]$data, paste(currentDir,outputFilename,sep = ""))

# Make UMAP by sample
f <- DimPlot(SCdata, reduction = "umap", group.by = "Sample", pt.size = 0.5, label = F) + 
  theme(text = element_text(size = 12)) 
show(f)
outputFilename <- "Control_UMAP_sample.csv"
write.csv(f[[1]]$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Make bar plot showing distribution of clusters
# ----------------------------------------------------------------

# Stacked bar plot showing clusters vs layer
f <- dittoBarPlot(SCdata, var = "seurat_clusters", group.by = "Sample")
show(f)
outputFilename <- "Control_cluster_sample.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Make bar plots for cell cycle distribution
# ----------------------------------------------------------------

# Stacked bar plot showing cell cycle vs cluster
f <- dittoBarPlot(SCdata, "cellcycle", group.by = "seurat_clusters")
show(f)
outputFilename <- "Control_cellcycle_cluster.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))

# Stacked bar plot showing cell cycle vs layer
f <- dittoBarPlot(SCdata, "cellcycle", group.by = "Sample")
show(f)
outputFilename <- "Control_cellcycle_sample.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Identify DEGs across (average expression) and explore pathways
# ----------------------------------------------------------------

# First get all DEGs
SCdata <- SetIdent(SCdata, value = "Sample")
all.markers_markers <- FindAllMarkers(SCdata, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

# Filter out any that have adjusted p-values > 0.05
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)

# Collect up unique DEGs
allDEGs <- unique(all.markers_markers$gene)
print(paste("There are",length(allDEGs),"DEGs between the different layers."))

# Calculate and plot average expression for all the DEGs
colormap = colorRamp2(c(-1.5, 0, 1.5), c("magenta", "black", "yellow"))
f <- dittoHeatmap(Avg_exp, genes = unique(all.markers_markers$gene),
                  name = "control-layer-DEGs",
                  scaled.to.max = FALSE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE, 
                  heatmap.colors = colormap,
                  fontsize = 7, complex = TRUE)
outputFilename <- "Control_allDEGs_sample.csv"
write.csv(f@matrix, paste(currentDir,outputFilename,sep = ""))
show(f)

# Get top 5 DEGs per cluster
all.markers_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
top5 <- top5 %>% as.data.frame() %>% arrange(cluster)

# Calculate and plot average expression for the top 5 DEGs
colormap = colorRamp2(c(-1.2, 0, 1.2), c("magenta", "black", "yellow"))
f <- dittoHeatmap(Avg_exp, genes = unique(top5$gene),
                  name = "control-layer-top5DEGs",
                  scaled.to.max = FALSE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE, 
                  heatmap.colors = colormap,
                  fontsize = 7, complex = TRUE)
outputFilename <- "Control_top5DEGs_sample.csv"
write.csv(f@matrix, paste(currentDir,outputFilename,sep = ""))
show(f)

# Clean up
rm(allDEGs, colormap)


# ----------------------------------------------------------------
# Explore pathways that differ between layers
# ----------------------------------------------------------------

# First get all DEGs
SCdata <- SetIdent(SCdata, value = "Sample")
all.markers_markers <- FindAllMarkers(SCdata, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

# Filter out any that have adjusted p-values > 0.05
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)
allDEGs <- unique(all.markers_markers$gene)

dfsample <- split(all.markers_markers$gene,all.markers_markers$cluster)

# Convert genes to ENTREZIDs
dfsample$`control-core` = bitr(dfsample$`control-core`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`control-middle` = bitr(dfsample$`control-middle`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`control-periphery` = bitr(dfsample$`control-periphery`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- list("control-core" = dfsample$`control-core`$ENTREZID,
                 "control-middle" = dfsample$`control-middle`$ENTREZID,
                 "control-periphery" = dfsample$`control-periphery`$ENTREZID)



### Alternatively, explore by cluster rather than layer

# First get all DEGs
SCdata <- SetIdent(SCdata, value = "seurat_clusters")
all.markers_markers <- FindAllMarkers(SCdata, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

# Filter out any that have adjusted p-values > 0.05
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)
allDEGs <- unique(all.markers_markers$gene)

dfsample <- split(all.markers_markers$gene,all.markers_markers$cluster)

# Convert genes to ENTREZIDs
dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- list("0" = dfsample$`0`$ENTREZID,
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID)


# Examine the most affected pathways
test_comparison <- compareCluster(genelist, fun = "enrichPathway")
f <- dotplot(test_comparison, x = "Cluster", color = "p.adjust", showCategory = 20,
        by = "geneRatio",
        size = "Count",
        split = NULL,
        includeAll = FALSE,
        font.size = 10,
        title = "",
        label_format = 50,
        group = FALSE,
        shape = FALSE)
outputFilename <- "Control_pathways_sample.csv"
write.csv(f[["data"]], paste(currentDir,outputFilename,sep = ""))
show(f)

# Have a look at gene ontology (biological processes)
test_comparison <- compareCluster(genelist, fun = "enrichGO", ont="BP", OrgDb = "org.Hs.eg.db")
f <- dotplot(test_comparison, x = "Cluster", color = "p.adjust", showCategory = 10,
             by = "geneRatio",
             size = "Count",
             split = NULL,
             includeAll = TRUE,
             font.size = 10,
             title = "",
             label_format = 50,
             group = FALSE,
             shape = FALSE)
show(f)

# Clean up
rm(all.markers_markers, test_comparison, dfsample, genelist)


# ----------------------------------------------------------------
# End of Figure 2, can delete things now
# ----------------------------------------------------------------

# Clear figure and output variables
rm(f, outputFilename)

# Clear SC experiment
rm(SCdata, Avg_exp)



# ================================================================
# Figures 4 and 5: Control and drug
# ================================================================


# ----------------------------------------------------------------
# Load the data
# ----------------------------------------------------------------

load("SCdata_drug_and_control_sbos.RData")


# ----------------------------------------------------------------
# Process the data
# ----------------------------------------------------------------

# Add some metadata about which samples received treatment
samples <- c(SBO03 = "Drug", SBO04 = "Vehicle",
             SBO05 = "Drug", SBO06 = "Vehicle",
             SBO07 = "Drug", SBO08 = "Vehicle")
names <- as.character(samples[SCdata$dscore])
SCdata <- AddMetaData(SCdata, metadata = names, col.name = 'Treatment')
rm(samples, names)

# Clustering
SCdata <- SetIdent(SCdata, value = SCdata$seurat_clusters)
SCdata <- FindClusters(SCdata, resolution = 0.4)
SCdata <- RunUMAP(SCdata, dims = 1:20)

# Renaming clusters
renamed_clusters <- c('0' = 4, '1' = 5, '2'= 1, '3' = 6, '4' = 2, '5' = 3, '6' = 9, '7' = 7, '8' = 8)
names <- as.factor(as.character(renamed_clusters[SCdata$seurat_clusters]))
SCdata <- AddMetaData(SCdata, metadata = names, col.name = 'renamed_clusters')

# Create a copy of the single-cell experiment for UMAP plots prior to subsetting
SCdata_all <- SCdata

# Create a separate dataset with cluster 6/9 removed
SCdata <- subset(x = SCdata, idents = 6, invert = TRUE)
normdata <- NormalizeData(SCdata, normalization.method = "LogNormalize", scale.factor = 10000)
normdata <- FindVariableFeatures(normdata, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(normdata)
scaled_data <- ScaleData(normdata)
dim_red <- RunPCA(scaled_data, features = VariableFeatures(object = scaled_data))
SCdata <- FindNeighbors(dim_red, dims = 1:20)

# Now do the clustering on this new object
SCdata <- FindClusters(SCdata, resolution = 0.4) # was 0.4, try smaller to reduce number of groups
SCdata <- RunUMAP(SCdata, dims = 1:20)

# Re-order cluster numbers; number in quotes is the original cluster, which is mapped to the corresponding integer
new_clusters <- c('0' = 4, '1' = 5, '2'= 1, '3' = 2, '4' = 3, '5' = 6, '6' = 7, '7' = 8, '8' = 9) 
SCdata$seurat_clusters <- as.factor(as.numeric(as.character(new_clusters[SCdata$seurat_clusters])))

# Clean up
rm(all.genes, names, new_clusters, renamed_clusters, normdata, scaled_data, dim_red)

# Scoring cell cycle
s.genes <- cc.genes.updated.2019$s.genes 
s.genes <- s.genes[!is.na(match(s.genes, rownames(SCdata)))]
g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- g2m.genes[!is.na(match(g2m.genes, rownames(SCdata)))] 
cellcycle = CellCycleScoring(SCdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) 
SCdata$cellcycle <- paste(cellcycle$Phase)

# Clean up
rm(s.genes, g2m.genes, cellcycle)

# Calculate average expression
Avg_exp_cluster <- AverageExpression(SCdata,
                                     assays = "RNA", return.seurat = TRUE,
                                     group.by = "seurat_clusters", slot = "data")
Avg_exp_sample <- AverageExpression(SCdata,
                                    assays = "RNA", return.seurat = TRUE,
                                    group.by = "Sample", slot = "data")


# ----------------------------------------------------------------
# Check metadata
# ----------------------------------------------------------------

# Subset by sample
SCdata <- SetIdent(SCdata, value = SCdata$Sample)
SCdata_controlcore <- subset(x = SCdata, idents = "control-core")
SCdata_controlmiddle <- subset(x = SCdata, idents = "control-middle")
SCdata_controlperiphery <- subset(x = SCdata, idents = "control-periphery")
SCdata_drugcore <- subset(x = SCdata, idents = "drug-core")
SCdata_drugmiddle <- subset(x = SCdata, idents = "drug-middle")
SCdata_drugperiphery <- subset(x = SCdata, idents = "drug-periphery")
listofsamples = c(SCdata_controlcore,SCdata_controlmiddle,SCdata_controlperiphery,
                  SCdata_drugcore,SCdata_drugmiddle,SCdata_drugperiphery,
                  SCdata)

# Metadata summary - transcript counts
for (sample in listofsamples) {
  print(summary(sample@meta.data[["nCount_RNA"]]))
}

# Metadata summary - gene counts
for (sample in listofsamples) {
  print(summary(sample@meta.data[["nFeature_RNA"]]))
}

# Metadata summary - SBO counts
for (sample in listofsamples) {
  print(summary(sample@meta.data[["nCount_HTO"]]))
}

# Clean up
rm(sample,listofsamples)
rm(SCdata_controlcore,SCdata_controlmiddle,SCdata_controlperiphery)
rm(SCdata_drugcore,SCdata_drugmiddle,SCdata_drugperiphery)


# ----------------------------------------------------------------
# Plot UMAPs
# ----------------------------------------------------------------

# According to clusters
f <- DimPlot(SCdata_all, reduction = "umap", group.by = "renamed_clusters",
             pt.size = 0.5, label = F) + theme(text = element_text(size = 12))
show(f)
outputFilename <- "All_UMAP_cluster.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))

# According to control/drug
f <- DimPlot(SCdata_all, reduction = "umap", group.by = "Treatment", cols = c("purple","green"),
             pt.size = 0.5, label = F) + theme(text = element_text(size = 12))
show(f)
outputFilename <- "All_UMAP_treatment.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Make bar plots for distributions
# ----------------------------------------------------------------

# Stacked bar plot showing cluster vs layer/treatment
f <- dittoBarPlot(SCdata, "renamed_clusters", group.by = "Sample")
show(f)
outputFilename <- "All_cluster_sample.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))

# Stacked bar plot showing cell cycle vs layer/treatment
f <- dittoBarPlot(SCdata, "cellcycle", group.by = "Sample")
show(f)
outputFilename <- "All_cellcycle_sample.csv"
write.csv(f$data, paste(currentDir,outputFilename,sep = ""))


# ----------------------------------------------------------------
# Single-cell differential gene expression
# ----------------------------------------------------------------

# Find all differentially expressed genes
SCdata <- SetIdent(SCdata, value = "seurat_clusters")
all.markers_markers <- FindAllMarkers(SCdata, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

# Filter out any that have adjusted p-values > 0.05 (doesn't really matter since we're only taking top 5)
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)

# Reorder as above and get top 5 per cluster
all.markers_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
top5 <- top5 %>% as.data.frame() %>% arrange(cluster)

# Plot heatmap
colormap = colorRamp2(c(-2, 0, 2), c("magenta", "black", "yellow"))
f <- dittoHeatmap(SCdata, genes = unique(top5$gene), annot.by = "seurat_clusters",
                  name = "all-layer-sc-DEGs",
                  scaled.to.max = FALSE,
                  use_raster = FALSE,
                  show_colnames = FALSE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE, 
                  heatmap.colors = colormap,
                  fontsize = 7, complex = TRUE)
f = draw(f)

# Save output
outputFilename <- "All_scheatmap_cluster.csv"
write.csv(f@ht_list[["all-layer-sc-DEGs"]]@matrix, paste(currentDir,outputFilename,sep = ""))

totals <- c()
for (i in 1:8) {
  totals <- c(totals, sum(SCdata$seurat_clusters == i))
}
outputFilename <- "All_scheatmap_cluster_totals.csv"
write.csv(totals, paste(currentDir,outputFilename,sep = ""))

# View by average expression as well for a quick look
f <- dittoHeatmap(Avg_exp_cluster, genes = unique(top5$gene),
                  name = "all-layer-sc-DEGs",
                  scaled.to.max = FALSE,
                  use_raster = FALSE,
                  show_colnames = FALSE,
                  cluster_cols = FALSE,
                  cluster_rows = FALSE, 
                  heatmap.colors = colormap,
                  fontsize = 7, complex = TRUE)
f = draw(f)

# Clean up
rm(i, totals, colormap, outputFilename)
rm(all.markers_markers, top5)


# ----------------------------------------------------------------
# Average differential gene expression
# ----------------------------------------------------------------

SCdata <- SetIdent(SCdata, value = SCdata$Sample)
all.markers_markers <- FindAllMarkers(SCdata, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

# Filter out any that have adjusted p-values > 0.05
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)

f <- dittoHeatmap(Avg_exp_sample, genes = all.markers_markers$gene,
                  name = "avg",
                  scaled.to.max = FALSE,
                  show_colnames = FALSE,
                  cluster_cols = FALSE,
                  cluster_rows = TRUE,
                  fontsize = 7, complex = TRUE)
f = draw(f)
show(f)

# Save output
outputFilename <- "Fig4f.csv"
write.csv(f@ht_list[["avg"]]@matrix, paste(currentDir,outputFilename,sep = ""))

# Also save order of rows
outputFilename <- "Fig4f_order.csv"
write.csv(row_order(f), paste(currentDir,outputFilename,sep = ""))

# Clean up
rm(all.markers_markers, f, colormap)


# ----------------------------------------------------------------
# Explore pathways that differ between layers (not really useful)
# ----------------------------------------------------------------

SCdata_all <- SetIdent(SCdata_all, value = "Sample")

# Get all DEGs if desired
drug.markers_markers <- FindAllMarkers(SCdata_all,
                                       min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)

# First get all DEGs in the drug samples
drug.markers_markers <- FindAllMarkers(
  subset(SCdata_all,
         ident = c("drug-core","drug-middle","drug-periphery")
         ), min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
drugDEGs <- unique(drug.markers_markers$gene)

# Now let's get DEGs in the control
control.markers_markers <- FindAllMarkers(subset(SCdata_all, ident = c("control-core","control-middle","control-periphery")),
                                          min.pct = 0.25, logfc.threshold = 0.1, return.thresh = 0.05)
control.markers_markers <- filter(control.markers_markers, p_val_adj < 0.05)
controlDEGs <- unique(control.markers_markers$gene)

# Subtract control from drug
drugonlyDEGs <- setdiff(drugDEGs, controlDEGs)
drug.markers_markers <- filter(drug.markers_markers, !gene %in% controlDEGs)

# Upregulated only
drug.markers_markers <- filter(drug.markers_markers, avg_log2FC > 0)
drug.markers_markers <- filter(drug.markers_markers, p_val_adj < 0.05)
drugDEGs <- unique(drug.markers_markers$gene)

# Average expression
Avg_exp <- AverageExpression(SCdata_all,
                             assays = "RNA",
                             return.seurat = TRUE,
                             group.by = "Sample",
                             slot = "data")
dfsample <- split(drug.markers_markers$gene,drug.markers_markers$cluster)

# Convert genes to ENTREZIDs
genelist <- c()
group_names <- names(dfsample)
for (x in group_names) {
  dfsample[[x]] = bitr(dfsample[[x]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genelist[[x]] <- dfsample[[x]]$ENTREZID
}
rm(x,dfsample,group_names)

# Examine the most affected pathways
#test_comparison <- compareCluster(genelist, fun = "enrichPathway")
test_comparison <- compareCluster(genelist, fun = "enrichGO", ont="BP", OrgDb = "org.Hs.eg.db")
f <- dotplot(test_comparison, x = "Cluster", showCategory = 8,
             color = "p.adjust", 
             by = "geneRatio",
             size = "Count",
             split = NULL,
             includeAll = TRUE,
             font.size = 8,
             title = "",
             label_format = 50,
             group = FALSE,
             shape = FALSE)
show(f)

# Clean up
rm(control.markers_markers, drug.markers_markers)
rm(controlDEGs, drugDEGs, allDEGs, drugonlyDEGs)
rm(test_comparison, genelist)



# ----------------------------------------------------------------
# Start Fig 5, spatially dependent genes
# ----------------------------------------------------------------

SCdata <- SetIdent(SCdata, value = "Sample")

# First get all DEGs in the drug samples and filter out any that have adjusted p-values > 0.05
drug.markers_markers <- FindAllMarkers(
  subset(SCdata,
         ident = c("drug-core","drug-middle","drug-periphery")
         ), min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
drug.markers_markers <- filter(drug.markers_markers, p_val_adj < 0.05)
drugDEGs <- unique(drug.markers_markers$gene)

# Now get DEGs in the control (these will be removed so we can set a loose threshold)
control.markers_markers <- FindAllMarkers(
  subset(SCdata,
         ident = c("control-core","control-middle","control-periphery")
         ),  min.pct = 0.25, logfc.threshold = 0.1, return.thresh = 0.05)
#control.markers_markers <- filter(control.markers_markers, p_val_adj < 0.05)
controlDEGs <- unique(control.markers_markers$gene)

# Subtract control from drug
drugonly.markers_markers <- filter(drug.markers_markers, !gene %in% controlDEGs)
drugonlyDEGs <- setdiff(drugDEGs, controlDEGs)

# Plot average expression
f <- dittoHeatmap(Avg_exp_sample,
                  genes = unique(drugonly.markers_markers$gene),
                  name = "avg",
                  scaled.to.max = FALSE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = TRUE, 
                  fontsize = 7, complex = TRUE)
f = draw(f)
show(f)

# Save output
outputFilename <- "Fig5a.csv"
write.csv(f@ht_list[["avg"]]@matrix, paste(currentDir,outputFilename,sep = ""))

# Clean up
rm(drug.markers_markers,control.markers_markers,drugonly.markers_markers)
rm(drugDEGs,controlDEGs,drugonlyDEGs)


# ----------------------------------------------------------------
# Spatially independent genes
# ----------------------------------------------------------------

# First get all DEGs between drug and untreated
SCdata <- SetIdent(SCdata, value = "Treatment")
treatment.markers_markers <- FindAllMarkers(
  subset(SCdata,
         ident = c("Drug","Vehicle")
  ), min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.05)
treatment.markers_markers <- filter(treatment.markers_markers, p_val_adj < 0.05)
treatmentDEGs <- unique(treatment.markers_markers$gene)

# Get DEGs between drug layers (these will be removed so we can set a loose threshold)
SCdata <- SetIdent(SCdata, value = "Sample")
drug.markers_markers <- FindAllMarkers(
  subset(SCdata,
         ident = c("drug-core","drug-middle","drug-periphery")
  ), min.pct = 0.25, logfc.threshold = 0.1, return.thresh = 0.05)
drug.markers_markers <- filter(drug.markers_markers, p_val_adj < 0.05)
drugDEGs <- unique(drug.markers_markers$gene)

# Now get DEGs in the control (these will be removed so we can set a loose threshold)
control.markers_markers <- FindAllMarkers(
  subset(SCdata,
         ident = c("control-core","control-middle","control-periphery")
  ),  min.pct = 0.25, logfc.threshold = 0.1, return.thresh = 0.05)
control.markers_markers <- filter(control.markers_markers, p_val_adj < 0.05)
controlDEGs <- unique(control.markers_markers$gene)

# Subtract drug and then control from treatment
nonspatial.markers_markers <- filter(treatment.markers_markers, !gene %in% drugDEGs)
nonspatial.markers_markers <- filter(nonspatial.markers_markers, !gene %in% controlDEGs)
nonspatial.markers_markers <- filter(nonspatial.markers_markers, avg_log2FC > 0)
nonspatial.markers_markers.control <- filter(nonspatial.markers_markers, cluster == "Vehicle")
nonspatial.markers_markers.drug <- filter(nonspatial.markers_markers, cluster == "Drug")
nonspatialDEGs <- setdiff(treatmentDEGs, drugDEGs)
nonspatialDEGs <- setdiff(nonspatialDEGs, controlDEGs)

# Plot average expression
f <- dittoHeatmap(Avg_exp_sample,
                  genes = unique(nonspatial.markers_markers.control$gene),
                  name = "Normalised\naverage\nexpression",
                  scaled.to.max = FALSE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = TRUE,
                  fontsize = 6.5, complex = TRUE)
f = draw(f)
show(f)
f <- dittoHeatmap(Avg_exp_sample,
                  genes = unique(nonspatial.markers_markers.drug$gene),
                  name = "Norm avg",
                  scaled.to.max = FALSE,
                  show_colnames = TRUE,
                  cluster_cols = FALSE,
                  cluster_rows = TRUE,
                  fontsize = 6.5, complex = TRUE)
f = draw(f)
show(f)

# Clean up
rm(treatment.markers_markers,drug.markers_markers,control.markers_markers,nonspatial.markers_markers)
rm(nonspatial.markers_markers.control,nonspatial.markers_markers.drug)
rm(drugDEGs,controlDEGs,nonspatialDEGs,treatmentDEGs)


# ----------------------------------------------------------------
# LINCS correlation analysis
# ----------------------------------------------------------------

# Want data organised by layer
Avg_exp <- Avg_exp_sample@assays[["RNA"]]@scale.data
Avg_exp <- cbind(symbol = rownames(Avg_exp), Avg_exp)
files <- c("LINCS/L1000_LINCS_DCIC_REP.A013_HELA_24H_D07_camptothecin_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A013_HELA_24H_D08_camptothecin_3.33uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A013_HELA_24H_D09_camptothecin_1.11uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.B013_HELA_24H_D09_camptothecin_0.25uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A013_HELA_24H_D11_camptothecin_0.125uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.B013_HELA_24H_D11_camptothecin_0.03uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.B013_HELA_24H_D12_camptothecin_0.01uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A005_HELA_24H_C13_irinotecan_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A005_HELA_24H_C14_irinotecan_3.33uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A005_HELA_24H_C15_irinotecan_1.11uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A005_HELA_24H_C16_irinotecan_0.37uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A026_HELA_24H_L19_doxorubicin_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A026_HELA_24H_L20_doxorubicin_3.33uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A026_HELA_24H_L21_doxorubicin_1.11uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.B026_HELA_24H_L21_doxorubicin_0.25uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.B026_HELA_24H_L24_doxorubicin_0.01uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A008_HELA_24H_A08_daunorubicin_3.33uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A008_HELA_24H_A09_daunorubicin_1.11uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A008_HELA_24H_A10_daunorubicin_0.37uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A008_HELA_24H_A11_daunorubicin_0.125uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A009_HELA_24H_H19_paclitaxel_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A011_HELA_24H_H13_docetaxel_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A009_HELA_24H_D19_chlorambucil_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A018_HELA_24H_P01_ifosfamide_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_REP.A004_HELA_24H_J19_fluorouracil_10uM.tsv",
           "LINCS/L1000_LINCS_DCIC_LJP008_HELA_24H_M13_decitabine_10uM.tsv")
samples = c("control-core","control-middle","control-periphery",
            "drug-core","drug-middle","drug-periphery")
corr <- c()
pvals <- c()

# Loop over each drug
for (file in files) {
  # Load file
  drug.signature <- read.delim(file, header = TRUE, sep = "\t", dec = ".")
  rownames(drug.signature) <- drug.signature$symbol
  # Merge dataframes so can correlate
  df_merge <- merge(drug.signature, Avg_exp, by="symbol")
  df_merge$symbol <- NULL
  # Calculate correlation with each layer
  corr_tau <- c()
  corr_p <- c()
  for (sample in samples)
  {
    corr_test <- cor.test(df_merge$CD.coefficient, as.numeric(df_merge[[sample]]), method = "kendall")
    corr_tau <- cbind(corr_tau, corr_test[["estimate"]][["tau"]])
    corr_p <- cbind(corr_p, corr_test[["p.value"]])
  }
  # Store result
  corr <- rbind(corr, corr_tau)
  pvals <- rbind(pvals, corr_p)
  print(paste("Completed file",nrow(corr),"of",length(files)))
}

# Clean up
rm(Avg_exp,files,file,drug.signature,df_merge,sample,corr_test,corr_tau,corr_p)

rownames(corr) <- c("10","3.33","1.11","0.25","0.125","0.03","0.01", # Camptothecin
                    "10","3.33","1.11","0.37", # Irinotecan
                    "10","3.33","1.11","0.25","0.01", # Doxorubicin
                    "3.33","1.11","0.37","0.125", # Daunorubicin
                    "Paclitaxel",
                    "Docetaxel",
                    "Chlorambucil",
                    "Ifosphamide",
                    "Fluorouracil",
                    "Decitabine"
)
colnames(corr) <- samples
rownames(pvals) <- rownames(corr)
colnames(pvals) <- samples

outputFilename <- "All_drug_correlation_r.csv"
write.csv(corr, paste(currentDir,outputFilename,sep = ""))
outputFilename <- "All_drug_correlation_p.csv"
write.csv(pvals, paste(currentDir,outputFilename,sep = ""))

col_fun = colorRamp2(c(-0.1, 0, 0.1), c("magenta", "black", "yellow"))
f <- Heatmap(corr, name = "drug_response_correlation",
             col = col_fun,
             cluster_rows = FALSE, cluster_columns = FALSE)
f = draw(f)

# Clean up
rm(corr,pvals,col_fun,samples)


# ----------------------------------------------------------------
# End of Figure 5, can delete things now
# ----------------------------------------------------------------

# Clear figure and output variables
rm(f, outputFilename)

# Clear SC experiment
rm(SCdata, SCdata_all, Avg_exp_cluster, Avg_exp_sample)













SCdata  <- SetIdent(SCdata, value = SCdata$Sample)
all.markers_markers <- FindAllMarkers(SCdata,
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25,
                                      return.thresh = 0.05)
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)

dfsample <- split(all.markers_markers$gene,all.markers_markers$cluster)
dfsample$`control-core` = bitr(dfsample$`control-core`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`control-middle` = bitr(dfsample$`control-middle`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`control-periphery` = bitr(dfsample$`control-periphery`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- list("control-core" = dfsample$`control-core`$ENTREZID,
                 "control-middle" = dfsample$`control-middle`$ENTREZID,
                 "control-periphery" = dfsample$`control-periphery`$ENTREZID)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
dotplot(GOclusterplot, x = "Cluster", color = "p.adjust", showCategory = 10, 
        by = "geneRatio",
        size = "Count",
        split = NULL,
        includeAll = TRUE,
        font.size = 10,
        title = "",
        label_format = 50,
        group = FALSE,
        shape = FALSE) 

# Alternatively:
all.markers_markers = FindMarkers(SCdata, ident.1 = "control-core", ident.2 = "control-periphery")
all.markers_markers = FindMarkers(SCdata, ident.1 = "drug-core", ident.2 = "drug-periphery")
all.markers_markers = FindMarkers(SCdata, ident.1 = "drug-periphery", ident.2 = "drug-core")
all.markers_markers <- filter(all.markers_markers, p_val_adj < 0.05)
dfsample <- data.frame(rownames(all.markers_markers), all.markers_markers$avg_log2FC, row.names = rownames(all.markers_markers))
colnames(dfsample) <- c("SYMBOL","avg_log2FC")
map <- bitr(rownames(dfsample), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample <- merge(map,dfsample)

geneList = dfsample[,3]
names(geneList) = as.character(dfsample[,2])
geneList = sort(geneList, decreasing = TRUE)

gse <- gseGO(geneList, ont = "ALL",
      OrgDb = "org.Hs.eg.db",
      keyType = "ENTREZID",
      exponent = 1,
      minGSSize = 1,
      maxGSSize = 500,
      eps = 1e-10,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = TRUE,
      seed = FALSE,
      by = "fgsea")
dotplot(gse,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 50,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        orderBy = "x",
        label_format = 100
)

