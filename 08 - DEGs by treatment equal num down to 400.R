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

mywd <- paste0("[path]/", myTime(), ' by treatment equalNum down to 400 cells')
dir.create(mywd)
setwd(mywd)

# Load downsampled data from "06 - Equalize_cell_number_down_to_400_cells.R"
load("[path]/s1p.combined_equalNum_down_to_400_cells.RData")

for (clusterArg in unique(s1p.combined.equalNum$lineageClusterIdents)) {
  patientList = c("s1p2", "s1p3", "s1p4", "s1p5", "s1p6", "s1p7", "s1p8")
  dayList = c("Day 1", "Day 3", "Day 7")
  
  DefaultAssay(s1p.combined.equalNum) <- "RNA"
  
  s1p.combined.equalNum <- NormalizeData(s1p.combined.equalNum)
  
  s1p.combined.equalNum <- SetIdent(s1p.combined.equalNum, value = "lineageClusterIdents")
  
  TestDEGsbyTreatment <- function(data, myClust, t) {
    myCluster <- subset(data, idents = myClust) %>% subset(subset = time == t)
    Idents(myCluster) <- "treatment"
    
    DEGs <- FindMarkers(myCluster, ident.1 = "BAF312", ident.2 = "Placebo", min.pct = 0.1)
  }
  FindDEGsbyTreatment <- function(data, myClust, t) {
  	myCluster <- subset(data, idents = myClust) %>% subset(subset = time == t)
  	Idents(myCluster) <- "treatment"
  	
  	DEGs <- FindMarkers(myCluster, ident.1 = "BAF312", ident.2 = "Placebo", logfc.threshold = 0, min.pct = 0.1)
  	
  	write.table(DEGs, file = paste0("DEGs by treatment equalNum down to 400 cells ", myClust, " ", t, ".txt"), 
  	            sep = "\t", col.names = NA)
  }
  
  print(paste("Testing for DEGs by treatment for", clusterArg))
  
  for (t in dayList) {
    possibleError <- tryCatch(TestDEGsbyTreatment(s1p.combined.equalNum, clusterArg, t),
                              error = function(e) {print(paste(i, t, "threw an error."))}
                              )
    if (!inherits(possibleError, "error")) {
      FindDEGsbyTreatment(s1p.combined.equalNum, clusterArg, t)
    }
  }
}
