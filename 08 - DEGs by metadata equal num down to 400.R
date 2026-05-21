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

mywd <- paste0("[path]/", myTime(), ' by metadata equalNum down to 400 cells')
dir.create(mywd)
setwd(mywd)

# Load downsampled data from "06 - Equalize_cell_number_down_to_400_cells.R"
load("[path]/s1p.combined_equalNum_down_to_400_cells.RData")
rm(s1p.combined.equalNum.BAF312)
rm(s1p.combined.equalNum.placebo)


# Calculate DEGs

patientList = c("s1p2", "s1p3", "s1p4", "s1p5", "s1p6", "s1p7", "s1p8")
dayList = c("Day 1", "Day 3", "Day 7")

DefaultAssay(s1p.combined.equalNum) <- "RNA"

s1p.combined.equalNum <- NormalizeData(s1p.combined.equalNum)

s1p.combined.equalNum <- SetIdent(s1p.combined.equalNum, value = "lineageClusterIdents")


TestDEGsbyMetadata <- function(data, myClust, t, metadataName, myContrast1, myContrast2) {
  myCluster <- subset(data, idents = myClust) %>% subset(subset = time == t)
  Idents(myCluster) <- metadataName
  
  DEGs <- FindMarkers(myCluster, ident.1 = myContrast1, ident.2 = myContrast2, min.pct = 0.1)
}
FindDEGsbyMetadata <- function(data, myClust, t, metadataName, myContrast1, myContrast2) {
  myCluster <- subset(data, idents = myClust) %>% subset(subset = time == t)
  Idents(myCluster) <- metadataName
  
  DEGs <- FindMarkers(myCluster, ident.1 = myContrast1, ident.2 = myContrast2, logfc.threshold = 0, min.pct = 0.1)
  
  write.table(DEGs, file = paste0("DEGs by metadata ", metadataName, " ", myContrast1, " - ", myContrast2, " equalNum down to 400 cells ", myClust, " ", t, ".txt"), 
              sep = "\t", col.names = NA)
}

# Build df of contrasts to test
contrastdf <- data.frame(contrastName = character(),
                         contrast1 = character(),
                         contrast2 = character())

contrastdf <- rbind(contrastdf, c('MRS_Day90_split_at_4', 'MRS_Day90_below_4', 'MRS_Day90_above_or_equal_to_4')) %>%
  rbind(c('NIHSS_Day1_split_at_10', 'NIHSS_Day1_above_or_equal_to_10', 'NIHSS_Day1_below_10')) %>%
  rbind(c('PHE_Day1_split_at_30', 'PHE_Day1_above_or_equal_to_30', 'PHE_Day1_below_30')) %>%
  rbind(c('Hematoma_Vol_Day1_split_at_30', 'Hematoma_Vol_Day1_above_or_equal_to_30', 'Hematoma_Vol_Day1_below_30')) %>%
  rbind(c('Rel_PHE_Day1_split_at_1', 'Rel_PHE_Day1_above_or_equal_to_1', 'Rel_PHE_Day1_below_1')) %>%
  rbind(c('ICH_score_split_at_2', 'ICH_score_above_or_equal_to_2', 'ICH_score_below_2'))

colnames(contrastdf) <- c('contrastName', 'contrast1', 'contrast2')
                      

for (clust in unique(s1p.combined.equalNum$lineageClusterIdents)) {
  writeLines(paste("\nTesting for DEGs by treatment for", clust))
  for (contrastNum in 1:nrow(contrastdf)) {
    contrastName <- contrastdf$contrastName[contrastNum]
    contrast1 <- contrastdf$contrast1[contrastNum]
    contrast2 <- contrastdf$contrast2[contrastNum]
    print(paste(contrastName, contrast1, contrast2))
    
    s1pTemp <- s1p.combined.equalNum
    
    for (t in dayList) {
        FindDEGsbyMetadata(s1pTemp, clust, t, contrastName, contrast1, contrast2)
    }
  }
}
