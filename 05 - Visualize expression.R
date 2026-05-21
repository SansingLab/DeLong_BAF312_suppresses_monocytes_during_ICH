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

mywd <- paste0("[path]/", myTime(), " sc expression total cells")
dir.create(mywd)
setwd(mywd)

# Load clustered data with metadata added from Script 03 - Add Patient Metadata
load("[path]/s1p.combined_metadata.RData")

s1p.combined <- NormalizeData(s1p.combined)

# # # Feature plots # # #
Idents(s1p.combined) <- 'lineageClusterIdents'
s1p.subsampled.10000 <- subset(s1p.combined, subset = (lineageClusterIdents != 'Dead' & 
                                                         lineageClusterIdents != 'Doublets'))[, sample(colnames(s1p.combined), size = 10000, replace=F)]
DefaultAssay(s1p.subsampled.10000) <- "RNA"

# Lineage markers for Fig 1
Idents(s1p.subsampled.10000) <- 'lineageClusterIdents'
FeaturePlot(s1p.subsampled.10000, pt.size = 1, ncol = 3, features = c("CD14", "FCGR3A", "NKG7", "MS4A1", "CD3E", "CD8A"))
ggsave(paste0("FeaturePlot lineage markers ", myTime(), ".pdf"), height = 7, width = 13)


# Dotplots for Supp Fig 1

DefaultAssay(s1p.combined) <- "RNA"

# Order broadClusters for DotPlot
s1p.combined$broadClusterIdents <- factor(s1p.combined$broadClusterIdents, 
                                        levels = rev(
                                          c('Monocytes',
                                            'Activated Monocytes',
                                            'CD16+ Monocytes',
                                            'cDC1', 'cDC2',
                                            'NK cells', 'NK CD56hi', 'NK prolif',
                                            'CD4 Naive', 'CD4 TCM', 'CD4 TEM',
                                            'Treg',
                                            'CD8 Naive', 'CD8 TCM', 'CD8 TEM',
                                            'MAIT', 'dnT',
                                            'B Naive', 'B Int', 'B Mem', 
                                            'Plasmablasts',
                                            'pDC',
                                            'HSPC',
                                            'Platelets',
                                            'Doublets',
                                            'Dead')
                                        ))
  
s1p.combined <- SetIdent(s1p.combined, value = 'broadClusterIdents')

pdf(file = paste0("DotPlot broadClusterIdents RNA Assay with Dead and Doublets ", myTime(), ".pdf"), height = 7, width = 12)
s1p.combined %>%
  DotPlot(features = c('ITGAM', 'CD14', # Monocytes
                       'IL1B', 'CCL3', # Activated Monos
                       'FCGR3A', # Non-classical monocytes
                       'HLA-DRB1', # MHCII
                       'CLEC9A', 'XCR1', # cDC1
                       'CD1C', 'CLEC10A', # cDC2
                       'GNLY', 'NKG7', 'GZMB', # NK cells
                       'KLRC1', # Naive (CD56hi) NK cells (Di Vito, 2019, Front. Imm.)
                       'MKI67',
                       'CD3E', # T cells 
                       'CCR7', # Naive T cells
                       'FOXP3', # Tregs
                       'CD8A',
                       'KLRB1', # NKT and MAIT cells
                       'MS4A1', 'IGHD', # B cells
                       'JCHAIN', 'CLEC4C', 'SERPINF1', # pDC
                       'RUNX1', 'HOXA9', # HPSC
                       'PPBP', #Platelets
                       'MT-CYB', 'MT-ATP6' # Dead cells
  ),
  cluster.idents = FALSE) + 
  RotatedAxis()
dev.off()
