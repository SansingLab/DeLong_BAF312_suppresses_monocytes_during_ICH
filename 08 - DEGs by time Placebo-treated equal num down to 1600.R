library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(sctransform)
library(EnhancedVolcano)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
        t <- gsub(":", "", Sys.time())
        t <- gsub(" ", "_", t)
        t <- gsub("\\.\\d+", "", t)
        return(t)
}

set.seed(1111)

mywd <- paste0("[path]/", myTime(), " Placebo by time equalNum down to 1600 cells")
dir.create(mywd)
setwd(mywd)

# Load downsampled data from "06 - Equalize_cell_number_down_to_1600_cells.R"
load("[path]/s1p.combined_equalNum_down_to_1600_cells.RData")

for (clusterArg in unique(s1p.combined.equalNum.Placebo$lineageClusterIdents)) {
  patientList = c("s1p2", "s1p3", "s1p4", "s1p5", "s1p6", "s1p7", "s1p8")
  dayList = c("Day 1", "Day 3", "Day 7")
  
  DefaultAssay(s1p.combined.equalNum.Placebo) <- "RNA"
  
  s1p.combined.equalNum.Placebo <- NormalizeData(s1p.combined.equalNum.Placebo)
  
  s1p.combined.equalNum.Placebo <- SetIdent(s1p.combined.equalNum.Placebo, value = "lineageClusterIdents")
  
  
  TestDEGsbyTime <- function(data, myClust) {
    myCluster <- subset(data, idents = myClust) %>% subset(subset = treatment == "Placebo")
    Idents(myCluster) <- "time"
    
    DEGs_D3vD1 <- FindMarkers(myCluster, ident.1 = "Day 3", ident.2 = "Day 1", logfc.threshold = 0, min.pct = 0.1)
    DEGs_D7vD1 <- FindMarkers(myCluster, ident.1 = "Day 7", ident.2 = "Day 1", logfc.threshold = 0, min.pct = 0.1)
    DEGs_D7vD3 <- FindMarkers(myCluster, ident.1 = "Day 7", ident.2 = "Day 3", logfc.threshold = 0, min.pct = 0.1)
  }
  FindDEGsbyTime <- function(data, myClust) {
    myCluster <- subset(data, idents = myClust) %>% subset(subset = treatment == "Placebo")
    Idents(myCluster) <- "time"
    
    # DEGs D3vD1
    DEGs_D3vD1 <- FindMarkers(myCluster, ident.1 = "Day 3", ident.2 = "Day 1", logfc.threshold = 0, min.pct = 0.1)
    write.table(DEGs_D3vD1, file = paste0("DEGs Placebo-treated equalNum down to 1600 cells ", myClust, " Day 3 vs Day 1.txt"), 
                sep = "\t", col.names = NA)
    
    # DEGs D7vD1
    DEGs_D7vD1 <- FindMarkers(myCluster, ident.1 = "Day 7", ident.2 = "Day 1", logfc.threshold = 0, min.pct = 0.1)
    write.table(DEGs_D7vD1, file = paste0("DEGs Placebo-treated equalNum down to 1600 cells ", myClust, " Day 7 vs Day 1.txt"), 
                sep = "\t", col.names = NA)
    
    # DEGs D7vD3
    DEGs_D7vD3 <- FindMarkers(myCluster, ident.1 = "Day 7", ident.2 = "Day 3", logfc.threshold = 0, min.pct = 0.1)
    write.table(DEGs_D7vD3, file = paste0("DEGs Placebo-treated equalNum down to 1600 cells ", myClust, " Day 7 vs Day 3.txt"), 
                sep = "\t", col.names = NA)
  }
  
  print(paste("Testing for DEGs by time for", clusterArg))
  
  possibleError <- tryCatch(TestDEGsbyTime(s1p.combined.equalNum.Placebo, clusterArg),
                            error = function(e) {print(paste0(clusterArg, "threw an error."))}
                            )
  if (!inherits(possibleError, 'error')) {
    FindDEGsbyTime(s1p.combined.equalNum.Placebo, clusterArg)
  }
}

