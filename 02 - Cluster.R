library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(sctransform)
library(limma)
# Libraries for SingleR
library(SingleR)
library(scuttle)

options(warn = 1) # set the warn level < 2 to prevent warnings from being converted into errors
options(bitmapType='cairo') # For showing graphics in RStudio Server

# Function to return current date and time in a filename-friendly format
myTime <- function() {
        t <- gsub(":", "", Sys.time())
        t <- gsub(" ", "_", t)
        t <- gsub("\\.\\d+", "", t)
        return(t)
}

mywd <- paste0("[path]/", myTime())
dir.create(mywd)
setwd(mywd)

set.seed(1111)
 
# Load integrated data from Script 01 - Integrate sctransform reference-based.R
load("[path]/integrate_sctransform_rpca_npcs100.RData")

# Scale the data to prepare for dimensional reduction
all.genes <- rownames(s1p.combined)

# Linear dimensional reduction (PCA)
s1p.combined <- ScaleData(s1p.combined, features = rownames(s1p.combined))
s1p.combined <- FindVariableFeatures(s1p.combined, selection.method = "vst")
s1p.combined <- RunPCA(s1p.combined, features = VariableFeatures(object = s1p.combined), npcs = 100)

Idents(s1p.combined) <- s1p.combined$time
DimPlot(s1p.combined, reduction = "umap")

DefaultAssay(s1p.combined) <- "RNA"

pdf(file = paste0("DimPlot Res 5 by lymphopenia SCT ", myTime(), ".pdf"), height = 7, width = 12)
DimPlot(s1p.combined, reduction = "umap", split.by = "lymphopenic", label = TRUE, raster = FALSE)
dev.off()


# Determine the dimensionality of the dataset
s1p.combined <- JackStraw(s1p.combined, num.replicate = 100, dims = 100)
gc()
s1p.combined <- ScoreJackStraw(s1p.combined, dims = 1:100)
#JackStrawPlot(s1p.combined, dims = 1:20)

pdf(file = paste0("ElbowPlot ", myTime(), ".pdf"))
ElbowPlot(s1p.combined, ndims = 100) 
dev.off()

pdf(file = paste0("Jackstraw JackStrawPlot 100 PCs ", myTime(), ".pdf"), height = 5, width = 13)
JackStrawPlot(object = s1p.combined, dims = 1:100)
dev.off()

pdf(file = paste0("Jackstraw JackStrawPlot 80 PCs ", myTime(), ".pdf"), height = 5, width = 12)
JackStrawPlot(object = s1p.combined, dims = 1:80)
dev.off()

save.image(file = "Clustering_JackStraw.RData")

# Cluster the cells
s1p.combined <- FindNeighbors(s1p.combined, dims = 1:50)
s1p.combined <- FindClusters(s1p.combined, resolution = 0.5) # Increase "resolution" to find more clusters. Authors recommend 0.6-1.2

# Non-linear dimentional reduction (UMAP)
s1p.combined <- RunUMAP(s1p.combined, dims = 1:50)

# Plot clusters
pdf(file = paste0("DimPlot Res 5 dims 50 ", myTime(), ".pdf"), height = 7, width = 10)
DimPlot(s1p.combined, reduction = "umap", label = TRUE)
dev.off()

pdf(file = paste0("DimPlot Res 5 dims 50 by lymphopenia ", myTime(), ".pdf"), height = 7, width = 12)
DimPlot(s1p.combined, reduction = "umap", split.by = "lymphopenic", label = TRUE)
dev.off()

pdf(file = paste0("DimPlot Res 5 dims 50 by patient ", myTime(), ".pdf"), height = 5, width =30)
DimPlot(s1p.combined, reduction = "umap", split.by = "patient", label = TRUE)
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
s1p.combined.markers <- FindAllMarkers(s1p.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers <- s1p.combined.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)

top10 <- s1p.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
s1p.combined.small <- subset(s1p.combined, downsample = 100)

png(file = paste0("Heatmap Res 5 ", myTime(), ".png"), res = 288, height = 7200, width = 5600)
DoHeatmap(s1p.combined.small, features = top10$gene)
dev.off()

save.image(file = "Res5_dims50_cluster_no_labels.RData")

# SingleR with Hao 2021 (Satija lab) dataset
# https://atlas.fredhutch.org/nygc/multimodal-pbmc/

DefaultAssay(s1p.combined) <- "RNA"

# Data from Hao 2021 downloaded from https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
hao <- LoadH5Seurat(file = "[path]/hao_multi.h5seurat") 

hao.subsampled <- hao[, sample(colnames(hao), size = 40000, replace = F)]
s1p.subsampled <- s1p.combined[, sample(colnames(s1p.combined), size = 40000, replace = F)]

s1psc <- as.SingleCellExperiment(s1p.subsampled) # Already less than 40,000 cells, don't need to subsample
haosc <- as.SingleCellExperiment(hao.subsampled)

s1pNorm <- logNormCounts(s1psc)
haoNorm <- logNormCounts(haosc)

remove(s1p.combined)
remove(s1psc)
remove(hao)
remove(haosc)
gc()

SingleR.results <- SingleR(test=s1pNorm, ref=haoNorm, labels=haoNorm$celltype.l2, de.method="wilcox")

table(SingleR.results$labels)

pdf(file = paste0("plotScoreHeatmap SingleR Hao Res 5 dims 50 ", myTime(), ".pdf"), height = 8, width =12)
plotScoreHeatmap(SingleR.results)
dev.off()

s1p.subsampled[["SingleR.labels"]] <- SingleR.results$labels

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels")
ggsave(paste0("./DimPlot SingleR Hao mapped onto S1P-combined with doublets res 5 ", myTime(), ".pdf"), width = 10, height = 7)

# Color only 8 clusters at a time for easier visualization
DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 0), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 1-8 mapped onto S1P-combined with doublets res 5 ", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 8), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 9-16 mapped onto S1P-combined with doublets res 5 ", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 16), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 17-24 mapped onto S1P-combined with doublets res 5 ", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 24), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 25-29 mapped onto S1P-combined with doublets res 5 ", myTime(), ".pdf"), width = 10, height = 7)
 

save.image(file = "Res5_dims50_cluster_no_labels_SingleR.RData")

# Label populations
load(file = "[path]/Res5_dims50_cluster_no_labels.RData")

new.cluster.ids <- c("Mono_1","Mono_2", "NK_1", "CD8_TEM_1", "Mono_CD16_1", "Mono_3",
                     "CD4_TCM_1", "CD4_Nv_1", "Mono_4", "B_Nv", "Mono-Lympho-Doublets",
                     "Mono_5", "cDC2", "NK_2", "B_Mem", "Mono_6",
                     "CD8_Nv_1", "CD4-CD8", "NK_CD56hi", "Treg", "Platelets",
                     "Plasmablast", "B_ Int", "MAIT", "pDC", "NK_prolif",
                     "HSPC", "cDC1", "CD8_TCM_2", "dnT")

names(new.cluster.ids) <- levels(s1p.combined.small)
s1p.combined.small <- RenameIdents(s1p.combined.small, new.cluster.ids) # Label the small subset
levels(s1p.combined.small) <- sort(new.cluster.ids) # Sort clusters in alphabetical order

png(file = paste0("Heatmap s1p-combined Res 5 Labels ", myTime(), ".png"), res = 288, height = 10000, width = 7000)
DoHeatmap(s1p.combined.small, features = top10$gene)
dev.off()

#names(new.cluster.ids) <- levels(s1p.combined)
s1p.combined <- RenameIdents(s1p.combined, new.cluster.ids) # Label the full dataset
#s1p.combined <- StashIdent(object = s1p.combined, save.name = "ClusterNames_0.6")

pdf(file = paste0("DimPlot Res 5 all ", myTime(), ".pdf"), height = 7, width = 10)
DimPlot(s1p.combined, reduction = "umap", label = TRUE)
dev.off()

pdf(file = paste0("DimPlot Res 5 by lymphopenia ", myTime(), ".pdf"), height = 7, width = 14)
DimPlot(s1p.combined, reduction = "umap", split.by = "lymphopenic", label = TRUE)
dev.off()

pdf(file = paste0("DimPlot Res 5 by time group by patient ", myTime(), ".pdf"), height = 7, width = 17)
DimPlot(s1p.combined, reduction = "umap", split.by = "time", label = TRUE, group.by = "patient")
dev.off()

pdf(file = paste0("DimPlot Res 5 by time ", myTime(), ".pdf"), height = 7, width = 17)
DimPlot(subset(s1p.combined, subset = patient == "s1p2"), split.by = "time", label = TRUE)
DimPlot(subset(s1p.combined, subset = patient == "s1p3"), split.by = "time", label = TRUE)
DimPlot(subset(s1p.combined, subset = patient == "s1p4"), split.by = "time", label = TRUE)
DimPlot(subset(s1p.combined, subset = patient == "s1p5"), split.by = "time", label = TRUE)
DimPlot(subset(s1p.combined, subset = patient == "s1p6"), split.by = "time", label = TRUE)
DimPlot(subset(s1p.combined, subset = patient == "s1p7"), split.by = "time", label = TRUE)
DimPlot(subset(s1p.combined, subset = patient == "s1p8"), split.by = "time", label = TRUE)
dev.off()

s1p.combined <- NormalizeData(s1p.combined, assay = "RNA")
Idents(s1p.combined) <- factor(Idents(s1p.combined), levels = rev(sort(levels(Idents(s1p.combined)))))

pdf(file = paste0("DotPlot s1p-combined Res 5 RNA Assay ", myTime(), ".pdf"), height = 8, width = 13)
DotPlot(s1p.combined, assay = "RNA",
        features = c('ITGAM', 'CD14', 'S100A8', # Monocytes
                     'FCGR3A', # Non-classical monocytes
                     'HLA-DRB1',
                     'CLEC9A', 'XCR1', 'C1orf54', # cDC1
                     'CD1C', 'FCER1A', 'CLEC10A', # DC2
                     'GNLY', 'NKG7', # NK cells
                     'KLRC1', # Naive (CD56hi) NK cells (Di Vito, 2019, Front. Imm.)
                     'KLRB1', # NKT cells
                     'CD3E', 'CD2', 'CD8A', # T cells 
                     'CCR7', # Naive T cells
                     'MS4A1', 'BANK1', 'IGHD', # B cells
                     'GZMB', 'JCHAIN', 'CLEC4C', # pDC
                     'PPBP', 'PF4' #Platelets
        ),
        cluster.idents = FALSE) + 
  RotatedAxis()
dev.off()


save.image(file = "s1p-combined_cluster_labels.RData")


# # # # # # # Skip the above clustering and labeling, load previous data and recluster at higher resolution # # # # # # # 

load(file = "[path]/s1p-combined_cluster_labels.RData")

# Recluster myeloid populations
s1p.myeloid <- subset(s1p.combined, idents = c("Mono_1","Mono_2", "Mono_CD16_1", "Mono_3", 
                                               "Mono_4", "Mono-Lympho-Doublets",
                                               "Mono_5", "cDC2", "Mono_6",
                                               "Platelets", "HSPC", "cDC1"))

# Linear dimensional reduction (PCA)
s1p.myeloid <- RunPCA(s1p.myeloid, features = VariableFeatures(object = s1p.myeloid), npcs = 100)

# Determine the dimensionality of the dataset
s1p.myeloid <- JackStraw(s1p.myeloid, num.replicate = 100, dims = 80)
s1p.myeloid <- ScoreJackStraw(s1p.myeloid, dims = 1:80)

pdf(file = paste0("ElbowPlot Myeloid 40 dims ", myTime(), ".pdf"))
ElbowPlot(s1p.myeloid, ndims = 40) # Elbow appears at PC10, suggesting most of the variance is contained in the first 10 PCs
dev.off()

pdf(file = paste0("JackStrawPlot Myeloid 40 PCs ", myTime(), ".pdf"), height = 5, width = 8)
JackStrawPlot(object = s1p.myeloid, dims = 1:40)
dev.off()

save(s1p.myeloid, file = "s1p-Myeloid_Unclustered.RData")
load(file = "[path]/s1p-Myeloid_Unclustered.RData")

# Cluster the cells
s1p.myeloid <- FindNeighbors(s1p.myeloid, dims = 1:20)
s1p.myeloid <- FindClusters(s1p.myeloid, resolution = 1.2) # Increase "resolution" to find more clusters. Authors recommend 0.6-1.2

# Non-linear dimentional reduction (UMAP/tSNE)
s1p.myeloid <- RunUMAP(s1p.myeloid, dims = 1:20)

pdf(file = paste0("DimPlot myeloid dims 20 Res 12 ", myTime(), ".pdf"), height = 7, width = 8)
DimPlot(s1p.myeloid, reduction = "umap", label = TRUE)
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
s1p.myeloid.markers <- FindAllMarkers(s1p.myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.myeloid.markers <- s1p.myeloid.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)

top10 <- s1p.myeloid.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

s1p.myeloid.small <- subset(s1p.myeloid, downsample = 100)

DefaultAssay(s1p.myeloid.small) <- "SCT"
png(file = paste0("Heatmap myeloid dims 20 Res 12 ", myTime(), ".png"), res = 288, 5000, width = 3500)
DoHeatmap(s1p.myeloid.small, features = top10$gene)
dev.off()

save(s1p.myeloid, 
     s1p.myeloid.small, 
     s1p.myeloid.markers, 
     cluster.myeloid.markers, 
     top10, 
     file = "s1p-Myeloid_Unlabeled_Res12.RData")
load(file = "[path]/s1p-Myeloid_Unlabeled_Res12.RData")


# # # Run SingleR to help annotate clusters # # #
DefaultAssay(s1p.myeloid) <- "RNA"


# SingleR with Mulder 2021 dataset for monocytes and DCs
mulder <- readRDS("[path]/2021_MNP_Verse.RDS")
mulder$mulderClusters <- mulder$Clusters

s1p.subsampled <- s1p.myeloid[, sample(colnames(s1p.myeloid), size = 40000, replace=F)]
mulder.NoMac <- subset(mulder,
                       mulderClusters == "CD16- Mono-16" |
                         mulderClusters == "CD16- Mono-21" |
                         mulderClusters == "CD16+ Mono" |
                         mulderClusters == "cDC1" |
                         mulderClusters == "DC2/DC3-12" |
                         mulderClusters == "DC2/DC3-15" |
                         mulderClusters == "mregDC" |
                         mulderClusters == "pre-DC" |
                         mulderClusters == "Proliferating cells")
mulder.subsampled <- mulder.NoMac[, sample(colnames(mulder.NoMac), size = 40000, replace=F)]

muldersc <- as.SingleCellExperiment(mulder.subsampled)
s1psc <- as.SingleCellExperiment(s1p.subsampled)

mulderNorm <- logNormCounts(muldersc)
s1pNorm <- logNormCounts(s1psc)

remove(muldersc)
remove(mulder)
gc()

SingleR.results <- SingleR(test=s1pNorm, ref=mulderNorm, labels=mulderNorm$Clusters, de.method="wilcox")
table(SingleR.results$labels)

s1p.subsampled[["SingleR.labels"]] <- SingleR.results$labels
DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels")
ggsave(paste0("./DimPlot SingleR Mulder mapped onto S1P-myeloid Res 12 ", myTime(), ".pdf"), width = 10, height = 7)

pdf(file = paste0("plotScoreHeatmap SingleR s1p-myeloid Mulder Res 12 dims 20 ", myTime(), ".pdf"), height = 8, width =12)
plotScoreHeatmap(SingleR.results)
dev.off()

save(SingleR.results,
     s1p.subsampled,
     file = "s1p-Myeloid_Res12_SingleR_Mulder.RData")


# SingleR with Hao 2021 (Satija lab) dataset
# This is the dataset used by Azimuth
# https://atlas.fredhutch.org/nygc/multimodal-pbmc/

# Data from Hao 2021 downloaded from https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
hao <- LoadH5Seurat(file = "[path]/hao_multi.h5seurat") 

hao.subsampled <- hao[, sample(colnames(hao), size = 40000, replace=F)]

s1psc <- as.SingleCellExperiment(s1p.subsampled)
haosc <- as.SingleCellExperiment(hao.subsampled)

s1pNorm <- logNormCounts(s1psc)
haoNorm <- logNormCounts(haosc)

remove(s1psc)
remove(hao)
remove(haosc)
gc()

SingleR.results <- SingleR(test=s1pNorm, ref=haoNorm, labels=haoNorm$celltype.l2, de.method="wilcox")

table(SingleR.results$labels)

s1p.subsampled[["SingleR.labels"]] <- SingleR.results$labels
DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels")
ggsave(paste0("./DimPlot SingleR Hao mapped onto S1P-myeloid", myTime(), ".pdf"), width = 10, height = 7)

# Color only 8 clusters at a time for easier visualization
DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 0), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 1-8 mapped onto S1P-myeloid", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 8), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 9-16 mapped onto S1P-myeloid", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 16), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 17-24 mapped onto S1P-myeloid", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 24), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 25-31 mapped onto S1P-myeloid", myTime(), ".pdf"), width = 10, height = 7)


pdf(file = paste0("plotScoreHeatmap SingleR s1p-myeloid Hao Res 12 dims 20 ", myTime(), ".pdf"), height = 8, width =12)
plotScoreHeatmap(SingleR.results)
dev.off()


save(SingleR.results,
     s1p.subsampled,
     file = "Myeloid_SingleR_Hao.RData")

# Label populations
new.cluster.ids.myeloid <- c("Mono_1", "Mono_2", "Mono_3", "Mono_4", "Mono_MHCIIhi_1", "Mono_5",
                             "Mono_6", "Mono_CD16_1", "Mono_CD16_2", "Mono_act", "Mono_CD16_3",
                             "Mono_IFN-resp", "Dead_Mono", "Mono_CD16_C1q", "Mono_7", "Doublet_Mono-T",
                             "Mono_8", "cDC2", "Mono_MHCIIhi_2", "Doublet", "Platelets",
                             "Dead_Mono_CD16", "HSPC", "cDC1", "Doublet_Mono_CD16-T", "Dead_cDC2")

s1p.myeloid.small <- subset(s1p.myeloid, downsample = 100)
names(new.cluster.ids.myeloid) <- levels(s1p.myeloid.small)
s1p.myeloid.small <- RenameIdents(s1p.myeloid.small, new.cluster.ids.myeloid) # Label the small subset
levels(s1p.myeloid.small) <- sort(new.cluster.ids.myeloid) # Sort clusters in alphabetical order

DefaultAssay(s1p.myeloid.small) <- "SCT"

png(file = paste0("./Heatmap myeloid Res 12 Labels ", myTime(), ".png"), res = 288, height = 7000, width = 5000)
DoHeatmap(s1p.myeloid.small, features = top10$gene)
dev.off()

pdf(file = paste0("DimPlot myeloid Res 12 Labels ", myTime(), ".pdf"), height = 7, width = 10)
DimPlot(s1p.myeloid.small, reduction = "umap", label = TRUE)
dev.off()

s1p.myeloid <- RenameIdents(s1p.myeloid, new.cluster.ids.myeloid) # Label the full dataset

pdf(file = paste0("DimPlot myeloid Res 12 Labels ", myTime(), ".pdf"), height = 7, width = 11)
DimPlot(s1p.myeloid, reduction = "umap", label = TRUE)
dev.off()

save(s1p.myeloid, file = "s1p-Myeloid_Labeled.RData")

remove(s1p.myeloid)
remove(s1p.myeloid.markers)
remove(s1p.myeloid.small)
remove(cluster.myeloid.markers)
remove(top10)

load(file = "[path]/s1p-Myeloid_Labeled.RData")


load(file = "[path]/s1p-combined_cluster_labels.RData")

# Recluster lymphoid populations, including doublet populations
s1p.lymphoid <- subset(s1p.combined, idents = c("NK_1", "CD8_TEM_1",
                                                "CD4_TCM_1", "CD4_Nv_1", "B_Nv",
                                                "NK_2", "B_Mem",
                                                "CD8_Nv_1", "CD4-CD8", "NK_CD56hi", "Treg",
                                                "Plasmablast", "B_ Int", "MAIT", "pDC", "NK_prolif",
                                                "CD8_TCM_2", "dnT"))

# Linear dimensional reduction (PCA)
s1p.lymphoid <- RunPCA(s1p.lymphoid, features = VariableFeatures(object = s1p.lymphoid), npcs = 80)

# Determine the dimensionality of the dataset
s1p.lymphoid <- JackStraw(s1p.lymphoid, num.replicate = 100, dims = 80)
s1p.lymphoid <- ScoreJackStraw(s1p.lymphoid, dims = 1:80)
#JackStrawPlot(s1p.lymphoid, dims = 1:20)

pdf(file = paste0("Elbow Plot lymphoid ", myTime(), ".pdf"))
ElbowPlot(s1p.lymphoid, ndims = 40) # Elbow appears at PC10, suggesting most of the variance is contained in the first 10 PCs
dev.off()

pdf(file = paste0("JackStrawPlot lymphoid 40 PCs ", myTime(), ".pdf"), height = 5, width = 8)
JackStrawPlot(object = s1p.lymphoid, dims = 1:40)
dev.off()

save(s1p.lymphoid, file = "s1p-Lymphoid_Unclustered.RData")
load("[path]/s1p-Lymphoid_Unclustered.RData")


# Cluster lymphoid cells at higher resolution
s1p.lymphoid <- FindNeighbors(s1p.lymphoid, dims = 1:40)
s1p.lymphoid <- FindClusters(s1p.lymphoid, resolution = 2.8) # Increase "resolution" to find more clusters. Authors recommend 0.6-1.2
#head(Idents(s1p.lymphoid), 5)

# Non-linear dimentional reduction (UMAP/tSNE)
s1p.lymphoid <- RunUMAP(s1p.lymphoid, dims = 1:40)

DimPlot(s1p.lymphoid, reduction = "umap", label = TRUE)
ggsave(paste0("./DimPlot lymphoid Res 28 ", myTime(), ".pdf"), width = 9, height = 7)

# Color only 8 clusters at a time for easier visualization
DimPlot(s1p.lymphoid, reduction = "umap", cols = c(rep("#999999", 0), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot lymphoid Res 28 clusters 0-7 ", myTime(), ".pdf"), width = 9, height = 7)

DimPlot(s1p.lymphoid, reduction = "umap", cols = c(rep("#999999", 8), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot lymphoid Res 28 clusters 8-15 ", myTime(), ".pdf"), width = 9, height = 7)

DimPlot(s1p.lymphoid, reduction = "umap", cols = c(rep("#999999", 16), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot lymphoid Res 28 clusters 16-23 ", myTime(), ".pdf"), width = 9, height = 7)

DimPlot(s1p.lymphoid, reduction = "umap", cols = c(rep("#999999", 24), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot lymphoid Res 28 clusters 24-31 ", myTime(), ".pdf"), width = 9, height = 7)

DimPlot(s1p.lymphoid, reduction = "umap", cols = c(rep("#999999", 32), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot lymphoid Res 28 clusters 32-39 ", myTime(), ".pdf"), width = 9, height = 7)

DimPlot(s1p.lymphoid, reduction = "umap", cols = c(rep("#999999", 40), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot lymphoid Res 28 clusters 40 ", myTime(), ".pdf"), width = 9, height = 7)

# find markers for every cluster compared to all remaining cells, report only the positive ones
s1p.lymphoid.markers <- FindAllMarkers(s1p.lymphoid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.lymphoid.markers <- s1p.lymphoid.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)

top10 <- s1p.lymphoid.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
s1p.lymphoid.small <- subset(s1p.lymphoid, downsample = 100)

png(file = paste0("Heatmap lymphoid Res 28 ", myTime(), ".png"), res = 288, height = 8000, width = 7000)
DoHeatmap(s1p.lymphoid.small, features = top10$gene)
dev.off()

save(s1p.lymphoid, 
     s1p.lymphoid.small, 
     s1p.lymphoid.markers, 
     cluster.lymphoid.markers, 
     top10, 
     file = "Lymphoid_Unlabeled_Res_28.RData")

load(file = "[path]/s1p-Lymphoid_Unlabeled_Res_28.RData")

unique(Idents(s1p.lymphoid))


# # # Run SingleR to help annotate clusters # # #

# SingleR with Hao 2021 (Satija lab) dataset
# This is the dataset used by Azimuth
# https://atlas.fredhutch.org/nygc/multimodal-pbmc/

DefaultAssay(s1p.lymphoid) <- "RNA"

# Data from Hao 2021 downloaded from https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
hao <- LoadH5Seurat(file = "[path]/hao_multi.h5seurat") 

s1p.subsampled <- s1p.lymphoid[, sample(colnames(s1p.lymphoid), size = 40000, replace=F)]
hao.subsampled <- hao[, sample(colnames(hao), size = 40000, replace=F)]

s1psc <- as.SingleCellExperiment(s1p.subsampled) # Already less than 40,000 cells, don't need to subsample
haosc <- as.SingleCellExperiment(hao.subsampled)

s1pNorm <- logNormCounts(s1psc)
haoNorm <- logNormCounts(haosc)

remove(s1psc)
remove(hao)
remove(haosc)
gc()

SingleR.results <- SingleR(test=s1pNorm, ref=haoNorm, labels=haoNorm$celltype.l2, de.method="wilcox")

table(SingleR.results$labels)

pdf(file = paste0("plotScoreHeatmap SingleR s1p-lymphoid Hao Res 12 dims 40 ", myTime(), ".pdf"), height = 8, width =12)
plotScoreHeatmap(SingleR.results)
dev.off()

s1p.subsampled[["SingleR.labels"]] <- SingleR.results$labels

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels")
ggsave(paste0("./DimPlot SingleR Hao mapped onto S1P-lymphoid res 12 ", myTime(), ".pdf"), width = 10, height = 7)

# Color only 8 clusters at a time for easier visualization
DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 0), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 1-8 mapped onto S1P-lymphoid res 12 ", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 8), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 9-16 mapped onto S1P-lymphoid res 12 ", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 16), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 17-24 mapped onto S1P-lymphoid res 12 ", myTime(), ".pdf"), width = 10, height = 7)

DimPlot(s1p.subsampled, reduction = "umap", group.by = "SingleR.labels", cols = c(rep("#999999", 24), brewer.pal(8, name = "Set1"), rep("#999999", 100)))
ggsave(paste0("./DimPlot SingleR Hao 25-29 mapped onto S1P-lymphoid res 12 ", myTime(), ".pdf"), width = 10, height = 7)

save(SingleR.results,
     s1p.subsampled,
     file = "s1p-Lymphoid_SingleR_Hao.RData")
 

# Label populations
new.cluster.ids.lymphoid <- c("NK_1", "CD8_TEM_1", "NK_2", "NK_3", "CD4_TCM_1", "CD4_Nv_1",
                              "CD4_Nv_2", "B_Nv_1", "B_Nv_2", "CD8_Nv", "NK_4",
                              "CD4_TCM_2", "CD4_Nv_3", "Dead_T_cell", "CD8_TEM_2", "CD8_TEM_3",
                              "CD8_TEM_4", "CD4_CTL", "CD4_TCM_3", "Dead_NK", "CD4_TEM",
                              "CD8_TCM", "Dead_CD8_TEM", "NK_5", "NK_CD56hi", "B_mem",
                              "Treg", "B_int_1", "Dead_B_Nv", "Plasmablast", "MAIT",
                              "pDC", "B_int_2", "NK_prolif", "Doublet_1", "CD4_TCM_4",
                              "Doublet_B-T", "dnT", "CD4_Nv_4", "Dead_pDC", "ILC")

# CD28 distinguishes CD4+ TCM from Tnv, per Shoiab Bukhari 2023 Cell Reports Medicine
# Low ribosomal protein transcripts in dead cells, per Ordonez-Rueda 2020 Cytometry
# TRAV1-2 and SLC4A10 are markers of CD8 MAITS, per Jijing Shi 2021 Front. Immunol.


DefaultAssay(s1p.lymphoid) <- "SCT"
DefaultAssay(s1p.lymphoid.small) <- "SCT"
VlnPlot(s1p.lymphoid, features = "FOXP3")

s1p.lymphoid.small <- subset(s1p.lymphoid, downsample = 100)
names(new.cluster.ids.lymphoid) <- levels(s1p.lymphoid.small)
s1p.lymphoid.small <- RenameIdents(s1p.lymphoid.small, new.cluster.ids.lymphoid) # Label the small subset
levels(s1p.lymphoid.small) <- sort(new.cluster.ids.lymphoid) # Sort clusters in alphabetical order

png(file = paste0("Heatmap lymphoid Res 28 Labels with doublets ", myTime(), ".png"), res = 288, height = 8000, width = 7000)
DoHeatmap(s1p.lymphoid.small, features = top10$gene)
dev.off()

pdf(file = paste0("DimPlot lymphoid Res 28 Labels with doublets ", myTime(), ".pdf"), height = 10, width = 12)
DimPlot(s1p.lymphoid.small, reduction = "umap", label = TRUE)
dev.off()


#names(new.cluster.ids) <- levels(s1p.combined)
s1p.lymphoid <- RenameIdents(s1p.lymphoid, new.cluster.ids.lymphoid) # Label the full dataset
#s1p.combined <- StashIdent(object = s1p.combined, save.name = "ClusterNames_0.6")

pdf(file = paste0("DimPlot lymphoid Res 28 Labels with doublets ", myTime(), ".pdf"), height = 10, width = 14)
DimPlot(s1p.lymphoid, reduction = "umap", label = TRUE)
dev.off()


save(s1p.myeloid, s1p.lymphoid, s1p.combined, file = 's1pObjectsLabeled.RData')

save(s1p.lymphoid, s1p.lymphoid.markers, cluster.lymphoid.markers, file = "s1pLymphoidLabeled.RData")

#save.image(file = "clustered_labeled.RData")
load(file = "[path]/s1pLymphoidLabeled.RData")


load(file = "[path]/s1pLymphoidLabeled.RData")
load(file = "[path]/s1p-Myeloid_Labeled.RData")
load(file = "[path]/s1p-combined_cluster_labels.RData")


# Adding subcluster idents to the combined object with the guidance from 
# https://bioinformatics.stackexchange.com/questions/18313/assigning-subcluster-idents-to-original-object

s1p.combined[["clusterIdents"]] <- as.character(Idents(s1p.combined))
s1p.lymphoid[["subclusterIdents"]] <- as.character(Idents(s1p.lymphoid))
s1p.myeloid[["subclusterIdents"]] <- as.character(Idents(s1p.myeloid))

new_metadata <- s1p.combined@meta.data
lymphoid_annotation <- s1p.lymphoid@meta.data
myeloid_annotation <- s1p.myeloid@meta.data

new_metadata$subclusterIdents <- sapply(rownames(new_metadata), 
                                        function(ita) ifelse(ita %in% rownames(lymphoid_annotation),
                                                             lymphoid_annotation[match(ita, rownames(lymphoid_annotation)), "subclusterIdents"],
                                                             new_metadata[match(ita, rownames(new_metadata)), "clusterIdents"]))

# For this second round of adding metadata, change last line to avoid writing over the first subclusterIdents
new_metadata$subclusterIdents <- sapply(rownames(new_metadata), 
                                        function(ita) ifelse(ita %in% rownames(myeloid_annotation),
                                                             myeloid_annotation[match(ita, rownames(myeloid_annotation)), "subclusterIdents"],
                                                             new_metadata[match(ita, rownames(new_metadata)), "subclusterIdents"]))

s1p.combined[["subclusterIdents"]] <- new_metadata$subclusterIdents

Idents(s1p.combined) <- s1p.combined@meta.data$subclusterIdents

# Order the subclusterIdents alphabetically
levels(s1p.combined) <- sort(unique(s1p.combined@meta.data$subclusterIdents))

DimPlot(s1p.combined, reduction = "umap", label = TRUE)

pdf(file = paste0("DimPlot subclustered with doublets Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 18)
DimPlot(s1p.combined, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()


# Create broader cluster idents to pool cell types
s1p.combined$broadClusterIdents <- NA
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "B_int_1" |
                                  s1p.combined$subclusterIdents == "B_int_2"] <- "B Int"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "B_mem"] <- "B Mem"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "B_Nv_1" |
                                  s1p.combined$subclusterIdents == "B_Nv_2"] <- "B Naive"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "CD4_CTL" |
                                  s1p.combined$subclusterIdents == "CD4_TEM"] <- "CD4 TEM"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "CD4_Nv_1" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_2" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_3" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_4"] <- "CD4 Naive"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "CD4_TCM_1" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_2" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_3" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_4"] <- "CD4 TCM"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "CD8_Nv"] <- "CD8 Naive"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "CD8_TCM"] <- "CD8 TCM"                                  
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "CD8_TEM_1" |
                                  s1p.combined$subclusterIdents == "CD8_TEM_2" |
                                  s1p.combined$subclusterIdents == "CD8_TEM_3" |
                                  s1p.combined$subclusterIdents == "CD8_TEM_4"] <- "CD8 TEM"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "cDC1"] <- "cDC1"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "cDC2"] <- "cDC2"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Dead_B_Nv" |
                                  s1p.combined$subclusterIdents == "Dead_CD8_TEM" |
                                  s1p.combined$subclusterIdents == "Dead_cDC2" |
                                  s1p.combined$subclusterIdents == "Dead_Mono" |
                                  s1p.combined$subclusterIdents == "Dead_Mono_CD16" |
                                  s1p.combined$subclusterIdents == "Dead_NK" |
                                  s1p.combined$subclusterIdents == "Dead_pDC" |
                                  s1p.combined$subclusterIdents == "Dead_T_cell"] <- "Dead"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "dnT"] <- "dnT"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Doublet" |
                                  s1p.combined$subclusterIdents == "Doublet_1" |
                                  s1p.combined$subclusterIdents == "Doublet_B-T" |
                                  s1p.combined$subclusterIdents == "Doublet_Mono_CD16-T" |
                                  s1p.combined$subclusterIdents == "Doublet_Mono-T"] <- "Doublets"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "HSPC"] <- "HSPC"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "MAIT"] <- "MAIT"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Mono_1" |
                                  s1p.combined$subclusterIdents == "Mono_2" |
                                  s1p.combined$subclusterIdents == "Mono_3" |
                                  s1p.combined$subclusterIdents == "Mono_4" |
                                  s1p.combined$subclusterIdents == "Mono_5" |
                                  s1p.combined$subclusterIdents == "Mono_6" |
                                  s1p.combined$subclusterIdents == "Mono_7" |
                                  s1p.combined$subclusterIdents == "Mono_8"] <- "Monocytes"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Mono_act" |
                                  s1p.combined$subclusterIdents == "Mono_IFN-resp" |
                                  s1p.combined$subclusterIdents == "Mono_MHCIIhi_1" |
                                  s1p.combined$subclusterIdents == "Mono_MHCIIhi_2"] <- "Activated Monocytes"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Mono_CD16_1" |
                                  s1p.combined$subclusterIdents == "Mono_CD16_2" |
                                  s1p.combined$subclusterIdents == "Mono_CD16_3" |
                                  s1p.combined$subclusterIdents == "Mono_CD16_C1q"] <- "CD16+ Monocytes"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "NK_1" |
                                  s1p.combined$subclusterIdents == "NK_2" |
                                  s1p.combined$subclusterIdents == "NK_3" |
                                  s1p.combined$subclusterIdents == "NK_4" |
                                  s1p.combined$subclusterIdents == "NK_5" |
                                  s1p.combined$subclusterIdents == "ILC"] <- "NK cells"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "NK_CD56hi"] <- "NK CD56hi"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "NK_prolif"] <- "NK prolif"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "pDC"] <- "pDC"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Plasmablast"] <- "Plasmablasts"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Platelets"] <- "Platelets"
s1p.combined$broadClusterIdents[s1p.combined$subclusterIdents == "Treg"] <- "Treg"

s1p.combined$broadClusterIdents <- factor(s1p.combined$broadClusterIdents, 
                                          levels = rev(
                                            c('MAIT',
                                              'Monocytes',
                                              'cDC1',
                                              'dnT',
                                              'NK cells', 
                                              'CD4 Naive',
                                              'B Naive', 
                                              'CD8 TEM',
                                              'Plasmablasts',
                                              'CD4 TEM',
                                              'Activated Monocytes',
                                              'NK CD56hi',
                                              'CD4 TCM',
                                              'Dead',
                                              'NK prolif',
                                              'CD8 TCM',
                                              'CD16+ Monocytes',
                                              'B Int', 
                                              'Treg',
                                              'Doublets',
                                              'B Mem', 
                                              'CD8 Naive',
                                              'pDC',
                                              'HSPC',
                                              'cDC2',
                                              'Platelets')
                                          ))

DefaultAssay(s1p.combined) <- "SCT"
VlnPlot(s1p.combined, features = "CD14", group.by = "subclusterIdents")
VlnPlot(s1p.combined, features = "SLC4A10", group.by = "subclusterIdents")
VlnPlot(s1p.combined, features = "CD8A", group.by = "subclusterIdents")

# Create lineage cluster idents to pool cell types
s1p.combined$lineageClusterIdents <- NA
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "B_int_1" |
                                  s1p.combined$subclusterIdents == "B_int_2" |
                                  s1p.combined$subclusterIdents == "B_mem" |
                                  s1p.combined$subclusterIdents == "B_Nv_1" |
                                  s1p.combined$subclusterIdents == "B_Nv_2"] <- "B cells"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "CD4_CTL" |
                                  s1p.combined$subclusterIdents == "CD4_TEM" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_1" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_2" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_3" |
                                  s1p.combined$subclusterIdents == "CD4_Nv_4" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_1" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_2" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_3" |
                                  s1p.combined$subclusterIdents == "CD4_TCM_4" |
                                  s1p.combined$subclusterIdents == "Treg"] <- "CD4+ T cells"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "CD8_Nv" |
                                  s1p.combined$subclusterIdents == "CD8_TCM" |                                  
                                  s1p.combined$subclusterIdents == "CD8_TEM_1" |
                                  s1p.combined$subclusterIdents == "CD8_TEM_2" |
                                  s1p.combined$subclusterIdents == "CD8_TEM_3" |
                                  s1p.combined$subclusterIdents == "CD8_TEM_4"] <- "CD8+ T cells"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "cDC1"] <- "cDC1"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "cDC2"] <- "cDC2"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "Dead_B_Nv" |
                                  s1p.combined$subclusterIdents == "Dead_CD8_TEM" |
                                  s1p.combined$subclusterIdents == "Dead_cDC2" |
                                  s1p.combined$subclusterIdents == "Dead_Mono" |
                                  s1p.combined$subclusterIdents == "Dead_Mono_CD16" |
                                  s1p.combined$subclusterIdents == "Dead_NK" |
                                  s1p.combined$subclusterIdents == "Dead_pDC" |
                                  s1p.combined$subclusterIdents == "Dead_T_cell"] <- "Dead"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "dnT"] <- "dnT cells"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "Doublet" |
                                  s1p.combined$subclusterIdents == "Doublet_1" |
                                  s1p.combined$subclusterIdents == "Doublet_B-T" |
                                  s1p.combined$subclusterIdents == "Doublet_Mono_CD16-T" |
                                  s1p.combined$subclusterIdents == "Doublet_Mono-T"] <- "Doublets"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "HSPC"] <- "HSPC"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "MAIT"] <- "MAIT"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "Mono_1" |
                                  s1p.combined$subclusterIdents == "Mono_2" |
                                  s1p.combined$subclusterIdents == "Mono_3" |
                                  s1p.combined$subclusterIdents == "Mono_4" |
                                  s1p.combined$subclusterIdents == "Mono_5" |
                                  s1p.combined$subclusterIdents == "Mono_6" |
                                  s1p.combined$subclusterIdents == "Mono_7" |
                                  s1p.combined$subclusterIdents == "Mono_8" |
                                  s1p.combined$subclusterIdents == "Mono_act" |
                                  s1p.combined$subclusterIdents == "Mono_IFN-resp" |
                                  s1p.combined$subclusterIdents == "Mono_MHCIIhi_1" |
                                  s1p.combined$subclusterIdents == "Mono_MHCIIhi_2"] <- "Monocytes"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "Mono_CD16_1" |
                                  s1p.combined$subclusterIdents == "Mono_CD16_2" |
                                  s1p.combined$subclusterIdents == "Mono_CD16_3" |
                                  s1p.combined$subclusterIdents == "Mono_CD16_C1q"] <- "CD16+ Monocytes"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "NK_1" |
                                  s1p.combined$subclusterIdents == "NK_2" |
                                  s1p.combined$subclusterIdents == "NK_3" |
                                  s1p.combined$subclusterIdents == "NK_4" |
                                  s1p.combined$subclusterIdents == "NK_5" |
                                  s1p.combined$subclusterIdents == "ILC" |
                                  s1p.combined$subclusterIdents == "NK_CD56hi" |
                                  s1p.combined$subclusterIdents == "NK_prolif" |
                                  s1p.combined$subclusterIdents == "ILC"] <- "NK cells"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "pDC"] <- "pDC"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "Plasmablast"] <- "Plasmablasts"
s1p.combined$lineageClusterIdents[s1p.combined$subclusterIdents == "Platelets"] <- "Platelets"

# Create broad lineage idents to distinguish between myeloid and lymphoid cells
s1p.combined$broadLineageIdents <- s1p.combined$lineageClusterIdents
s1p.combined$broadLineageIdents[s1p.combined$lineageClusterIdents == "pDC" |
                                  s1p.combined$lineageClusterIdents == "Plasmablasts" |
                                  s1p.combined$lineageClusterIdents == "B cells" |
                                  s1p.combined$lineageClusterIdents == "dnT cells" |
                                  s1p.combined$lineageClusterIdents == "MAIT" |
                                  s1p.combined$lineageClusterIdents == "CD8+ T cells" |
                                  s1p.combined$lineageClusterIdents == "CD4+ T cells" |
                                  s1p.combined$lineageClusterIdents == "NK cells"] <- "Lymphoid"
s1p.combined$broadLineageIdents[s1p.combined$lineageClusterIdents == "cDC2" |
                                  s1p.combined$lineageClusterIdents == "cDC1" |
                                  s1p.combined$lineageClusterIdents == "cDC2" |
                                  s1p.combined$lineageClusterIdents == "CD16+ Monocytes" |
                                  s1p.combined$lineageClusterIdents == "Monocytes"] <- "Myeloid"

pdf(file = paste0("DimPlot lineageClusterIdents no doublets or dead Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 15)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead")) %>%
  DimPlot(reduction = "umap", label = TRUE, group.by = "lineageClusterIdents", raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot lineageClusterIdents no doublets or dead No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 15)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead")) %>%
  DimPlot(reduction = "umap", label = FALSE, group.by = "lineageClusterIdents", raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot broadClusterIdents no doublets or dead Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 16)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead")) %>%
  DimPlot(reduction = 'umap', label = TRUE, group.by = 'broadClusterIdents', raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot broadClusterIdents no doublets or dead No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 16)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead")) %>%
  DimPlot(reduction = 'umap', label = FALSE, group.by = 'broadClusterIdents', raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot subclusterIdents no doublets or dead Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 17)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead")) %>%
  DimPlot(reduction = 'umap', label = TRUE, group.by = 'subclusterIdents', raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot subclusterIdents no doublets or dead No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 17)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead")) %>%
  DimPlot(reduction = 'umap', label = FALSE, group.by = 'subclusterIdents', raster = FALSE)
dev.off()


pdf(file = paste0("DimPlot lineageClusterIdents Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 15)
  DimPlot(s1p.combined, reduction = "umap", label = TRUE, group.by = "lineageClusterIdents", raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot lineageClusterIdents No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 15)
  DimPlot(s1p.combined, reduction = "umap", label = FALSE, group.by = "lineageClusterIdents", raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot broadClusterIdents Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 16)
  DimPlot(s1p.combined, reduction = 'umap', label = TRUE, group.by = 'broadClusterIdents', raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot broadClusterIdents No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 16)
  DimPlot(s1p.combined, reduction = 'umap', label = FALSE, group.by = 'broadClusterIdents', raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot subclusterIdents Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 17)
  DimPlot(s1p.combined, reduction = 'umap', label = TRUE, group.by = 'subclusterIdents', raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot subclusterIdents No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 17)
  DimPlot(s1p.combined, reduction = 'umap', label = FALSE, group.by = 'subclusterIdents', raster = FALSE)
dev.off()


pdf(file = paste0("DimPlot broadClusterIdents Dead & Doublets ONLY No Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 16)
DimPlot(s1p.combined, reduction = "umap", label = FALSE, group.by = "broadClusterIdents", raster = FALSE, 
        cols = c(rep("#999999", 6), "#DD3333", rep("#999999", 5), "#3333DD", rep("#999999", 100)))
dev.off()

pdf(file = paste0("DimPlot broadClusterIdents Dead & Doublets ONLY Labels non-rasterized ", myTime(), ".pdf"), height = 10, width = 16)
DimPlot(s1p.combined, reduction = 'umap', label = TRUE, group.by = 'broadClusterIdents', raster = FALSE,
        cols = c(rep("#999999", 6), "#DD3333", rep("#999999", 5), "#3333DD", rep("#999999", 100)))
dev.off()


# DimPlots by patient
pdf(file = paste0("DimPlot lineageClusterIdents by patient Day 3 Placebo non-rasterized ", myTime(), ".pdf"), height = 6, width = 12)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead" & treatment == "Placebo" & time == "Day 3")) %>%
  DimPlot(reduction = "umap", label = TRUE, group.by = "lineageClusterIdents", split.by = "patient", raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot lineageClusterIdents by patient Day 3 BAF312 non-rasterized ", myTime(), ".pdf"), height = 6, width = 20)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead" & treatment == "BAF312" & time == "Day 3")) %>%
  DimPlot(reduction = "umap", label = TRUE, group.by = "lineageClusterIdents", split.by = "patient", raster = FALSE)
dev.off()

pdf(file = paste0("DimPlot lineageClusterIdents Day 3 BAF312 non-rasterized ", myTime(), ".pdf"), height = 6, width = 32)
subset(s1p.combined, subset = (broadClusterIdents != 'Doublets' & broadClusterIdents != "Dead" & time == "Day 1")) %>%
  DimPlot(reduction = "umap", label = TRUE, group.by = "lineageClusterIdents", split.by = "patient", raster = FALSE)
dev.off()

save(s1p.combined, file = "s1p.combined_subclustered.RData")

load("[path]/s1p.combined_subclustered.RData")# Color only 8 clusters at a time for easier visualization

load("[path]/s1p.combined_subclustered.RData")

