library(Seurat)
library(cowplot)
library(tidyverse)
library(RColorBrewer)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]", myTime(), " by time")
dir.create(mywd)
setwd(mywd)

DotPlotSelectedGenesPlacebo <- function(cluster, myFeatures, geneset = "select genes") {
  mySubset <- subset(s1p.combined.equalNum, treatment == 'Placebo' & lineageClusterIdents == cluster)
  
  (DotPlot(mySubset, features = myFeatures,
           group.by = "time",
           scale = TRUE) +
      scale_size(breaks = c(0, 20, 40, 60, 80, 100)) +
      coord_flip())
  ggsave(filename = paste0('Dotplot by time placebo ', cluster, ' scDEGs ', geneset, ' scaled ', myTime(), ".pdf"), width = 5, height = 1 + 0.2 * length(myFeatures))
}

DotPlotSelectedGenesBAF312 <- function(cluster, myFeatures, geneset = "select genes") {
  mySubset <- subset(s1p.combined.equalNum, treatment == 'BAF312' & lineageClusterIdents == cluster)
  
  (DotPlot(mySubset, features = myFeatures,
           group.by = "time",
           scale = TRUE) +
      scale_size(breaks = c(0, 20, 40, 60, 80, 100)) +
      coord_flip())
  ggsave(filename = paste0('Dotplot by time BAF312 ', cluster, ' scDEGs ', geneset, ' scaled ', myTime(), ".pdf"), width = 5, height = 1 + 0.2 * length(myFeatures))
}

# # # # # # # # # # # # # # # # # # # # # #
# Downsampled to 1,600 cells (monocytes)  #
# # # # # # # # # # # # # # # # # # # # # #

# Load downsampled cells from script 06 - Equalize cell number down to 1600 cells.R
load("[path]/s1p.combined_equalNum_down_to_1600_cells.RData")
rm(s1p.combined.equalNum.BAF312)
rm(s1p.combined.equalNum.placebo)
gc()

s1p.combined.equalNum <- NormalizeData(s1p.combined.equalNum)

DefaultAssay(s1p.combined.equalNum) <- "RNA"

myFeatures <- rev(c(
  'TNFAIP2',
  'NFKBIA',
  'FOSB',
  'JUN',
  'CD83',
  'EGR1',
  'CXCL8',
  'IRF1',
  'IRF7',
  'IL4R',
  'ITGA5',
  'VEGFA',
  'CLU', 
  'LGALS3', 
  'TGFB1',
  'ABCA1',
  'LDLR',
  'SREBF2'
))

DotPlotSelectedGenesPlacebo('Monocytes', myFeatures)
DotPlotSelectedGenesBAF312('Monocytes', myFeatures)



