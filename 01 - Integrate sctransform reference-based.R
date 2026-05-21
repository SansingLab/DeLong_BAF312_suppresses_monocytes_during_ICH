library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(sctransform)

options(warn = 1) # set the warn level < 2 to prevent warnings from being converted into errors

myTime <- function() {
	t <- gsub(":", "", Sys.time())
	t <- gsub(" ", "_", t)
	t <- gsub("\\.\\d+", "", t)
	return(t)
}

mywd <- paste0("[path]/", myTime())
dir.create(mywd)
setwd(mywd)

# Load aligned 10X data
s1p2d1.data <- Read10X(data.dir = "[path]/JHD1922_D1_Total_filtered_feature_bc_matrix")
s1p2d3.data <- Read10X(data.dir = "[path]/JHD1922_D3_Total_filtered_feature_bc_matrix")
s1p2d7.data <- Read10X(data.dir = "[path]/JHD1922_D7_Total_filtered_feature_bc_matrix")
s1p5d1.data <- Read10X(data.dir = "[path]/JHD1929_D1_filtered_feature_bc_matrix")
s1p5d3.data <- Read10X(data.dir = "[path]/JHD1929_D3_filtered_feature_bc_matrix")
s1p5d7.data <- Read10X(data.dir = "[path]/JHD1929_D7_filtered_feature_bc_matrix")
s1p34mixA.data <- Read10X(data.dir = "[path]/JHD1941_MIXA_filtered_feature_bc_matrix")
s1p34mixB.data <- Read10X(data.dir = "[path]/JHD1941_MIXB_filtered_feature_bc_matrix")
s1p34mixC.data <- Read10X(data.dir = "[path]/JHD1941_MIXC_filtered_feature_bc_matrix")
s1p78mixA.data <- Read10X(data.dir = "[path]/JHD2021_MIXA_filtered_feature_bc_matrix")
s1p78mixB.data <- Read10X(data.dir = "[path]/JHD2021_MIXB_filtered_feature_bc_matrix")
s1p78mixC.data <- Read10X(data.dir = "[path]/JHD2021_MIXC_filtered_feature_bc_matrix")
s1p6mixA.data <- Read10X(data.dir = "[path]/JHD2025_MIXA_filtered_feature_bc_matrix")
s1p6mixB.data <- Read10X(data.dir = "[path]/JHD2025_MIXB_filtered_feature_bc_matrix")

s1p2d1 <- CreateSeuratObject(counts = s1p2d1.data, min.cells = 3, min.features = 200)
s1p2d1$time <- "Day 1"
s1p2d1$patient <- "s1p2"
s1p2d1$lymphopenic <- "No"
s1p2d1$treatment <- "BAF312_interrupted"

s1p2d3 <- CreateSeuratObject(counts = s1p2d3.data, min.cells = 3, min.features = 200)
s1p2d3$time <- "Day 3"
s1p2d3$patient <- "s1p2"
s1p2d3$lymphopenic <- "No"
s1p2d3$treatment <- "BAF312_interrupted"

s1p2d7 <- CreateSeuratObject(counts = s1p2d7.data, min.cells = 3, min.features = 200)
s1p2d7$time <- "Day 7"
s1p2d7$patient <- "s1p2"
s1p2d7$lymphopenic <- "No"
s1p2d3$treatment <- "BAF312_interrupted"


s1p5d1 <- CreateSeuratObject(counts = s1p5d1.data, min.cells = 3, min.features = 200)
s1p5d1$time <- "Day 1"
s1p5d1$patient <- "s1p5"
s1p5d1$lymphopenic <- "No"
s1p5d1$treatment <- "Placebo"

s1p5d3 <- CreateSeuratObject(counts = s1p5d3.data, min.cells = 3, min.features = 200)
s1p5d3$time <- "Day 3"
s1p5d3$patient <- "s1p5"
s1p5d3$lymphopenic <- "No"
s1p5d3$treatment <- "Placebo"

s1p5d7 <- CreateSeuratObject(counts = s1p5d7.data, min.cells = 3, min.features = 200)
s1p5d7$time <- "Day 7"
s1p5d7$patient <- "s1p5"
s1p5d7$lymphopenic <- "No"
s1p5d7$treatment <- "Placebo"

s1p34mixA <- CreateSeuratObject(counts = s1p34mixA.data$'Gene Expression')
s1p34mixA[['HTO']] = CreateAssayObject(counts = s1p34mixA.data$'Antibody Capture')
s1p34mixB <- CreateSeuratObject(counts = s1p34mixB.data$'Gene Expression')
s1p34mixB[['HTO']] = CreateAssayObject(counts = s1p34mixB.data$'Antibody Capture')
s1p34mixC <- CreateSeuratObject(counts = s1p34mixC.data$'Gene Expression')
s1p34mixC[['HTO']] = CreateAssayObject(counts = s1p34mixC.data$'Antibody Capture')

s1p78mixA <- CreateSeuratObject(counts = s1p78mixA.data$'Gene Expression')
s1p78mixA[['HTO']] = CreateAssayObject(counts = s1p78mixA.data$'Antibody Capture')
s1p78mixB <- CreateSeuratObject(counts = s1p78mixB.data$'Gene Expression')
s1p78mixB[['HTO']] = CreateAssayObject(counts = s1p78mixB.data$'Antibody Capture')
s1p78mixC <- CreateSeuratObject(counts = s1p78mixC.data$'Gene Expression')
s1p78mixC[['HTO']] = CreateAssayObject(counts = s1p78mixC.data$'Antibody Capture')

s1p6mixA <- CreateSeuratObject(counts = s1p6mixA.data$'Gene Expression')
s1p6mixA[['HTO']] = CreateAssayObject(counts = s1p6mixA.data$'Antibody Capture')
s1p6mixB <- CreateSeuratObject(counts = s1p6mixB.data$'Gene Expression')
s1p6mixB[['HTO']] = CreateAssayObject(counts = s1p6mixB.data$'Antibody Capture')

# Remove unnecessary objects from environment
remove(s1p2d1.data)
remove(s1p2d3.data)
remove(s1p2d7.data)
remove(s1p5d1.data)
remove(s1p5d3.data)
remove(s1p5d7.data)
remove(s1p34mixA.data)
remove(s1p34mixB.data)
remove(s1p34mixC.data)
remove(s1p78mixA.data)
remove(s1p78mixB.data)
remove(s1p78mixC.data)
remove(s1p6mixA.data)
remove(s1p6mixB.data)

# Normalize the hashed data: s1p3, s1p4
s1p34mixA <- NormalizeData(s1p34mixA, assay = 'HTO', normalization.method = 'CLR')
s1p34mixB <- NormalizeData(s1p34mixB, assay = 'HTO', normalization.method = 'CLR')
s1p34mixC <- NormalizeData(s1p34mixC, assay = 'HTO', normalization.method = 'CLR')

# Demultiplex data: s1p3, s1p4
s1p34mixA <- HTODemux(s1p34mixA, assay = 'HTO', positive.quantile = 0.99)
s1p34mixB <- HTODemux(s1p34mixB, assay = 'HTO', positive.quantile = 0.99)
s1p34mixC <- HTODemux(s1p34mixC, assay = 'HTO', positive.quantile = 0.99)

# Sort HTO Idents to make the ridge plots more organized
Idents(s1p34mixA) <- factor(Idents(s1p34mixA), levels = rev(sort(levels(Idents(s1p34mixA)))))
Idents(s1p34mixB) <- factor(Idents(s1p34mixB), levels = rev(sort(levels(Idents(s1p34mixB)))))
Idents(s1p34mixC) <- factor(Idents(s1p34mixC), levels = rev(sort(levels(Idents(s1p34mixC)))))
 

# Visualize the results of demultiplexing
RidgePlot(s1p34mixA, assay = 'HTO',
          features = c('anti-human-Hashtag-1', 'anti-human-Hashtag-2', 
                       'anti-human-Hashtag-3', 'anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p34mixA ", myTime(), ".pdf"), height = 10, width = 10)

RidgePlot(s1p34mixB, assay = 'HTO',
          features = c('anti-human-Hashtag-1', 'anti-human-Hashtag-2', 
                       'anti-human-Hashtag-3', 'anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p34mixB ", myTime(), ".pdf"), height = 10, width = 10)

RidgePlot(s1p34mixC, assay = 'HTO',
          features = c('anti-human-Hashtag-1', 'anti-human-Hashtag-2', 
                       'anti-human-Hashtag-3', 'anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p34mixC ", myTime(), ".pdf"), height = 10, width = 10)

# Check HTO tags
table(s1p34mixA$HTO_classification.global)
table(s1p34mixB$HTO_classification.global)
table(s1p34mixC$HTO_classification.global)

# Check for doublets visually
FeatureScatter(s1p34mixA, feature1 = "anti-human-Hashtag-1", feature2 = "anti-human-Hashtag-2")

# Group cells based on the max HTO signal (Not sure how this works)
Idents(s1p34mixA) <- "HTO_maxID"
Idents(s1p34mixA) <- factor(Idents(s1p34mixA), levels = rev(sort(levels(Idents(s1p34mixA)))))
Idents(s1p34mixB) <- "HTO_maxID"
Idents(s1p34mixB) <- factor(Idents(s1p34mixB), levels = rev(sort(levels(Idents(s1p34mixB)))))
Idents(s1p34mixC) <- "HTO_maxID"
Idents(s1p34mixC) <- factor(Idents(s1p34mixC), levels = rev(sort(levels(Idents(s1p34mixC)))))

# Allows you to subset on HTO features
Idents(s1p34mixA) <- "HTO_classification.global"
Idents(s1p34mixB) <- "HTO_classification.global"
Idents(s1p34mixC) <- "HTO_classification.global"

VlnPlot(s1p34mixA, features = 'nCount_RNA', pt.size = 0.1, log = TRUE)

FeatureScatter(s1p34mixA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Remove cells that stained with no HTOs (negative cells)
s1p34mixA <- subset(s1p34mixA, idents = "Negative", invert = TRUE)
s1p34mixB <- subset(s1p34mixB, idents = "Negative", invert = TRUE)
s1p34mixC <- subset(s1p34mixC, idents = "Negative", invert = TRUE)

# Calculate t-SNE, plot
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p34mixA, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p34mixA <- RunTSNE(s1p34mixA, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p34mixA, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p34mixA')
ggsave(filename = paste0("Demultiplex Dimplot s1p34mixA ", myTime(), ".pdf"), height = 4, width = 5)

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p34mixB, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p34mixB <- RunTSNE(s1p34mixB, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p34mixB, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p34mixB')
ggsave(filename = paste0("Demultiplex Dimplot s1p34mixB ", myTime(), ".pdf"), height = 4, width = 5)

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p34mixC, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p34mixC <- RunTSNE(s1p34mixC, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p34mixC, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p34mixC')
ggsave(filename = paste0("Demultiplex Dimplot s1p34mixC ", myTime(), ".pdf"), height = 4, width = 5)

print("# # # # # Finished plotting t-SNE of s1p34mixA # # # # #")

# Subset on singlets
s1p34mixA <- subset(s1p34mixA, idents = "Singlet")
s1p34mixB <- subset(s1p34mixB, idents = "Singlet")
s1p34mixC <- subset(s1p34mixC, idents = "Singlet")

print("# # # # # Finished subsetting on singlets # # # # #")

# Merge the hashed samples together (merge instead of integrate to avoid repeated normalization)
# merge.data specifies whether to include normalized data
s1p34.combined <- merge(s1p34mixA, y = c(s1p34mixB, s1p34mixC), merge.data = TRUE) 

print("# # # # # Finished merging the hashed samples together # # # # #")

# Split the hashed samples into their original samples
Idents(s1p34.combined) <- "HTO_maxID"

s1p3d1 <- subset(s1p34.combined, idents = "anti-human-Hashtag-1")
s1p3d1$time <- "Day 1"
s1p3d1$patient <- "s1p3"
s1p3d1$lymphopenic <- "Yes"
s1p3d1$treatment <- "BAF312"

s1p3d3 <- subset(s1p34.combined, idents = "anti-human-Hashtag-2")
s1p3d3$time <- "Day 3"
s1p3d3$patient <- "s1p3"
s1p3d3$lymphopenic <- "Yes"
s1p3d3$treatment <- "BAF312"

s1p3d7 <- subset(s1p34.combined, idents = "anti-human-Hashtag-3")
s1p3d7$time <- "Day 7"
s1p3d7$patient <- "s1p3"
s1p3d7$lymphopenic <- "Yes"
s1p3d7$treatment <- "BAF312"

s1p4d1 <- subset(s1p34.combined, idents = "anti-human-Hashtag-4")
s1p4d1$time <- "Day 1"
s1p4d1$patient <- "s1p4"
s1p4d1$lymphopenic <- "Yes"
s1p4d1$treatment <- "BAF312"

s1p4d3 <- subset(s1p34.combined, idents = "anti-human-Hashtag-5")
s1p4d3$time <- "Day 3"
s1p4d3$patient <- "s1p4"
s1p4d3$lymphopenic <- "Yes"
s1p4d3$treatment <- "BAF312"

s1p4d7 <- subset(s1p34.combined, idents = "anti-human-Hashtag-6")
s1p4d7$time <- "Day 7"
s1p4d7$patient <- "s1p4"
s1p4d7$lymphopenic <- "Yes"
s1p4d7$treatment <- "BAF312"

print("# # # # # Finished splitting the hashed samples back into their original samples # # # # #")

# Normalize the hashed data: s1p7, s1p8
s1p78mixA <- NormalizeData(s1p78mixA, assay = 'HTO', normalization.method = 'CLR')
s1p78mixB <- NormalizeData(s1p78mixB, assay = 'HTO', normalization.method = 'CLR')
s1p78mixC <- NormalizeData(s1p78mixC, assay = 'HTO', normalization.method = 'CLR')

# Demultiplex data: s1p7, s1p8
s1p78mixA <- HTODemux(s1p78mixA, assay = 'HTO', positive.quantile = 0.99)
s1p78mixB <- HTODemux(s1p78mixB, assay = 'HTO', positive.quantile = 0.99)
s1p78mixC <- HTODemux(s1p78mixC, assay = 'HTO', positive.quantile = 0.99)

# Sort HTO Idents to make the ridge plots more organized
Idents(s1p78mixA) <- factor(Idents(s1p78mixA), levels = rev(sort(levels(Idents(s1p78mixA)))))
Idents(s1p78mixB) <- factor(Idents(s1p78mixB), levels = rev(sort(levels(Idents(s1p78mixB)))))
Idents(s1p78mixC) <- factor(Idents(s1p78mixC), levels = rev(sort(levels(Idents(s1p78mixC)))))

# Visualize the results of demultiplexing
RidgePlot(s1p78mixA, assay = 'HTO', 
          features = c('anti-human-Hashtag-1', 'anti-human-Hashtag-2', 
                       'anti-human-Hashtag-3', 'anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p78mixA ", myTime(), ".pdf"), height = 10, width = 10)

RidgePlot(s1p78mixB, assay = 'HTO', 
          features = c('anti-human-Hashtag-1', 'anti-human-Hashtag-2', 
                       'anti-human-Hashtag-3', 'anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p78mixB ", myTime(), ".pdf"), height = 10, width = 10)

RidgePlot(s1p78mixC, assay = 'HTO', 
          features = c('anti-human-Hashtag-1', 'anti-human-Hashtag-2', 
                       'anti-human-Hashtag-3', 'anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p78mixC ", myTime(), ".pdf"), height = 10, width = 10)

# Check HTO tags
table(s1p78mixA$HTO_classification.global)
table(s1p78mixB$HTO_classification.global)
table(s1p78mixC$HTO_classification.global)

# Group cells based on the max HTO signal (Not sure how this works)
Idents(s1p78mixA) <- "HTO_maxID"
Idents(s1p78mixB) <- "HTO_maxID"
Idents(s1p78mixC) <- "HTO_maxID"

# Check for doublets visually
FeatureScatter(s1p78mixA, feature1 = "anti-human-Hashtag-1", feature2 = "anti-human-Hashtag-2")

# Allows you to subset on HTO features
Idents(s1p78mixA) <- "HTO_classification.global"
Idents(s1p78mixB) <- "HTO_classification.global"
Idents(s1p78mixC) <- "HTO_classification.global"

VlnPlot(s1p78mixA, features = 'nCount_RNA', pt.size = 0.1, log = TRUE)

FeatureScatter(s1p78mixA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Remove cells that stained with no HTOs (negative cells)
s1p78mixA <- subset(s1p78mixA, idents = "Negative", invert = TRUE)
s1p78mixB <- subset(s1p78mixB, idents = "Negative", invert = TRUE)
s1p78mixC <- subset(s1p78mixC, idents = "Negative", invert = TRUE)

# Calculate t-SNE, plot
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p78mixA, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p78mixA <- RunTSNE(s1p78mixA, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p78mixA, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p78mixA')
ggsave(filename = paste0("Demultiplex Dimplot s1p78mixA ", myTime(), ".pdf"), height = 4, width = 5) 

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p78mixB, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p78mixB <- RunTSNE(s1p78mixB, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p78mixB, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p78mixB')
ggsave(filename = paste0("Demultiplex Dimplot s1p78mixB ", myTime(), ".pdf"), height = 4, width = 5) 

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p78mixC, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p78mixC <- RunTSNE(s1p78mixC, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p78mixC, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p78mixC')
ggsave(filename = paste0("Demultiplex Dimplot s1p78mixC ", myTime(), ".pdf"), height = 4, width = 5) 


print("# # # # # Finished plotting t-SNE of s1p78mixA # # # # #")

# Subset on singlets
Idents(s1p78mixA) <- "HTO_classification.global"
Idents(s1p78mixB) <- "HTO_classification.global"
Idents(s1p78mixC) <- "HTO_classification.global"
s1p78mixA <- subset(s1p78mixA, idents = "Singlet")
s1p78mixB <- subset(s1p78mixB, idents = "Singlet")
s1p78mixC <- subset(s1p78mixC, idents = "Singlet")

print("# # # # # Finished subsetting on singlets # # # # #")

# Merge the hashed samples together (merge instead of integrate to avoid repeated normalization)
# merge.data specifies whether to include normalized data
s1p78.combined <- merge(s1p78mixA, y = c(s1p78mixB, s1p78mixC), merge.data = TRUE) 

print("# # # # # Finished merging the hashed samples together # # # # #")

# Split the hashed samples into their original samples
Idents(s1p78.combined) <- "HTO_maxID"

s1p7d1 <- subset(s1p78.combined, idents = "anti-human-Hashtag-1")
s1p7d1$time <- "Day 1"
s1p7d1$patient <- "s1p7"
s1p7d1$lymphopenic <- "Yes"
s1p7d1$treatment <- "BAF312"

s1p7d3 <- subset(s1p78.combined, idents = "anti-human-Hashtag-2")
s1p7d3$time <- "Day 3"
s1p7d3$patient <- "s1p7"
s1p7d3$lymphopenic <- "Yes"
s1p7d3$treatment <- "BAF312"

s1p7d7 <- subset(s1p78.combined, idents = "anti-human-Hashtag-3")
s1p7d7$time <- "Day 7"
s1p7d7$patient <- "s1p7"
s1p7d7$lymphopenic <- "Yes"
s1p7d7$treatment <- "BAF312"

s1p8d1 <- subset(s1p78.combined, idents = "anti-human-Hashtag-4")
s1p8d1$time <- "Day 1"
s1p8d1$patient <- "s1p8"
s1p8d1$lymphopenic <- "Yes"
s1p8d1$treatment <- "BAF312"

s1p8d3 <- subset(s1p78.combined, idents = "anti-human-Hashtag-5")
s1p8d3$time <- "Day 3"
s1p8d3$patient <- "s1p8"
s1p8d3$lymphopenic <- "Yes"
s1p8d3$treatment <- "BAF312"

s1p8d7 <- subset(s1p78.combined, idents = "anti-human-Hashtag-6")
s1p8d7$time <- "Day 7"
s1p8d7$patient <- "s1p8"
s1p8d7$lymphopenic <- "Yes"
s1p8d7$treatment <- "BAF312"

print("# # # # # Finished splitting the hashed samples back into their original samples # # # # #")

# Normalize the hashed data: s1p6
s1p6mixA <- NormalizeData(s1p6mixA, assay = 'HTO', normalization.method = 'CLR')
s1p6mixB <- NormalizeData(s1p6mixB, assay = 'HTO', normalization.method = 'CLR')

# Remove hashtags 1-3, which were included in cellranger by mistake
#s1p6mixA@assays$HTO@counts <- s1p6mixA@assays$HTO@counts[4:6,]

# Demultiplex data: s1p6
s1p6mixA <- HTODemux(s1p6mixA, assay = 'HTO', positive.quantile = 0.99)
s1p6mixB <- HTODemux(s1p6mixB, assay = 'HTO', positive.quantile = 0.99)

# Sort HTO Idents to make the ridge plots more organized
Idents(s1p6mixA) <- factor(Idents(s1p6mixA), levels = rev(sort(levels(Idents(s1p6mixA)))))
Idents(s1p6mixB) <- factor(Idents(s1p6mixB), levels = rev(sort(levels(Idents(s1p6mixB)))))

# Visualize the results of demultiplexing
RidgePlot(s1p6mixA, assay = 'HTO', 
          features = c('anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p6mixA ", myTime(), ".pdf"), height = 6, width = 10)

RidgePlot(s1p6mixB, assay = 'HTO', 
          features = c('anti-human-Hashtag-4', 
                       'anti-human-Hashtag-5', 'anti-human-Hashtag-6'),
          ncol = 2)
ggsave(filename = paste0("HTO RidgePlot s1p6mixB ", myTime(), ".pdf"), height = 6, width = 10)

# Check HTO tags
table(s1p6mixA$HTO_classification.global)
table(s1p6mixB$HTO_classification.global)

# Group cells based on the max HTO signal (Not sure how this works)
Idents(s1p6mixA) <- "HTO_maxID"
Idents(s1p6mixB) <- "HTO_maxID"

# Check for doublets visually
FeatureScatter(s1p6mixA, feature1 = "anti-human-Hashtag-4", feature2 = "anti-human-Hashtag-5")

# Allows you to subset on HTO features
Idents(s1p6mixA) <- "HTO_classification.global"
Idents(s1p6mixB) <- "HTO_classification.global"

VlnPlot(s1p6mixA, features = 'nCount_RNA', pt.size = 0.1, log = TRUE)

FeatureScatter(s1p6mixA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Remove cells that stained with no HTOs (negative cells)
s1p6mixA <- subset(s1p6mixA, idents = "Negative", invert = TRUE)
s1p6mixB <- subset(s1p6mixB, idents = "Negative", invert = TRUE)

# Calculate t-SNE, plot
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p6mixA, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p6mixA <- RunTSNE(s1p6mixA, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p6mixA, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p6mixA')
ggsave(filename = paste0("Demultiplex Dimplot s1p6mixA ", myTime(), ".pdf"), height = 4, width = 5) 

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = s1p6mixB, assay = 'HTO')))) # Calculate a distance matrix using HTO
s1p6mixB <- RunTSNE(s1p6mixB, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(s1p6mixB, cols = c("Singlet" = "midnightblue", "Doublet" = "maroon1")) +
  labs(title = 's1p6mixB')
ggsave(filename = paste0("Demultiplex Dimplot s1p6mixB ", myTime(), ".pdf"), height = 4, width = 5) 

print("# # # # # Finished plotting t-SNE of s1p6mixA # # # # #")

# Subset on singlets
Idents(s1p6mixA) <- "HTO_classification.global"
Idents(s1p6mixB) <- "HTO_classification.global"
s1p6mixA <- subset(s1p6mixA, idents = "Singlet")
s1p6mixB <- subset(s1p6mixB, idents = "Singlet")

print("# # # # # Finished subsetting on singlets # # # # #")

# Merge the hashed samples together (merge instead of integrate to avoid repeated normalization)
# merge.data specifies whether to include normalized data
s1p6.combined <- merge(s1p6mixA, y = s1p6mixB, merge.data = TRUE) 

print("# # # # # Finished merging the hashed samples together # # # # #")

# Split the hashed samples into their original samples
Idents(s1p6.combined) <- "HTO_maxID"

s1p6d1 <- subset(s1p6.combined, idents = "anti-human-Hashtag-4")
s1p6d1$time <- "Day 1"
s1p6d1$patient <- "s1p6"
s1p6d1$lymphopenic <- "No"
s1p6d1$treatment <- "Placebo"

s1p6d3 <- subset(s1p6.combined, idents = "anti-human-Hashtag-5")
s1p6d3$time <- "Day 3"
s1p6d3$patient <- "s1p6"
s1p6d3$lymphopenic <- "No"
s1p6d3$treatment <- "Placebo"

s1p6d7 <- subset(s1p6.combined, idents = "anti-human-Hashtag-6")
s1p6d7$time <- "Day 7"
s1p6d7$patient <- "s1p6"
s1p6d7$lymphopenic <- "No"
s1p6d7$treatment <- "Placebo"

print("# # # # # Finished splitting the hashed samples back into their original samples # # # # #")

# Make a list of all the Seurat objects to make manipulation easier
s1p.list <- list(s1p2d1, s1p2d3, s1p2d7, 
                 s1p3d1, s1p3d3, s1p3d7, 
                 s1p4d1, s1p4d3, s1p4d7, 
                 s1p5d1, s1p5d3, s1p5d7,
                 s1p6d1, s1p6d3, s1p6d7,
                 s1p7d1, s1p7d3, s1p7d7,
                 s1p8d1, s1p8d3, s1p8d7)

# Remove mixed objects from environment
remove(s1p34mixA)
remove(s1p34mixB)
remove(s1p34mixC)
remove(s1p78mixA)
remove(s1p78mixB)
remove(s1p78mixC)
remove(s1p6mixA)
remove(s1p6mixB)
remove(s1p34.combined)
remove(s1p6.combined)
remove(s1p78.combined)
remove(hto.dist.mtx)

for (i in 1:length(s1p.list)) {
  # Set default assay to RNA
  DefaultAssay(s1p.list[[i]]) <- "RNA"
  
  # Add % MT gene stat
  s1p.list[[i]][["percent.mt"]] <- PercentageFeatureSet(s1p.list[[i]], pattern = "^MT-")
  
}

print("# # # # # Finished adding MT stat # # # # #")

# Visualize QC metrics as a violin plot
for (i in s1p.list) {
  Idents(i) <- paste(i$patient, i$time)
  VlnPlot(i, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("QC metrics ", unique(i$patient), " ", unique(i$time), " ", myTime(), ".pdf"), height = 4, width = 6)
}

# Remove individual objects from environment
remove(s1p2d1, s1p2d3, s1p2d7, 
       s1p3d1, s1p3d3, s1p3d7, 
       s1p4d1, s1p4d3, s1p4d7, 
       s1p5d1, s1p5d3, s1p5d7,
       s1p6d1, s1p6d3, s1p6d7,
       s1p7d1, s1p7d3, s1p7d7,
       s1p8d1, s1p8d3, s1p8d7)

# Filter dead cells
for (i in 1:length(s1p.list)) {
  s1p.list[[i]] <- subset(s1p.list[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25)
}

# Visualize QC metrics as a violin plot after filtering out dead cells
for (i in s1p.list) {
  #Idents(i) <- paste(i$patient, i$time)
  VlnPlot(i, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("QC metrics after filtering ", unique(i$patient), " ", unique(i$time), " ", myTime(), ".pdf"), height = 4, width = 6)
}

# Normalize with SCTransform
for (i in 1:length(s1p.list)) {
  s1p.list[[i]] <- SCTransform(s1p.list[[i]], vars.to.regress = "percent.mt")
}

print("# # # # # Finished normalizing # # # # #")

DefaultAssay(s1p.combined) <- "integrated"

# Integrate the datasets together with rpca
for (i in 1:length(s1p.list)) {
  s1p.list[[i]] <- FindVariableFeatures((s1p.list[[i]]))
}

s1p.features <- SelectIntegrationFeatures(object.list = s1p.list)

s1p.list <- lapply(X = s1p.list, FUN = function(x) {
  x <- ScaleData(x, features = s1p.features, verbose = FALSE)
  x <- RunPCA(x, features = s1p.features, verbose = FALSE, npcs = 100)
})

s1p.anchors <- FindIntegrationAnchors(object.list = s1p.list, reduction = "rpca", dims = 1:100)

s1p.combined <- IntegrateData(anchorset = s1p.anchors, dims = 1:100)

s1p.combined <- ScaleData(s1p.combined, verbose = FALSE)
s1p.combined <- RunPCA(s1p.combined, verbose = FALSE, npcs = 100)
s1p.combined <- RunUMAP(s1p.combined, dims = 1:50)

save.image(file = "integrate_sctransform_rpca_npcs100.RData")





