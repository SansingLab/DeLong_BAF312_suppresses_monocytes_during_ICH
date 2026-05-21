library(Seurat)
library(tidyverse)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]/", myTime(), " equalized cell number down to 400 cells setseed")
dir.create(mywd)
setwd(mywd)

# Load clustered data with metadata added from Script 03 - Add Patient Metadata
load("[path]/s1p.combined_metadata.RData")

myTable <- function(x) {
  return(table(x$lineageClusterIdents, x$patient, x$time))
}

subsample <- function(x) {
  normIDs <- vector()
  for (c in unique(x$lineageClusterIdents)) {
    print(c)
    s1p.clust <- subset(x, lineageClusterIdents == c)

    for (p in unique(s1p.clust$patient_time_lineage)) {
      ncells <- ncol(s1p.clust[, s1p.clust$patient_time_lineage == p])
      print(paste(p, ncells))
      if (ncells < 400) {
        normIDs <- c(normIDs, colnames(s1p.clust[, s1p.clust$patient_time_lineage == p]))
        print(paste('ncells for', p, ncells))
      } else {
        myCols <- colnames(s1p.clust[, s1p.clust$patient_time_lineage == p])
        set.seed(2222)
        newNormIDs <- sample(myCols, size = 400, replace = FALSE)
        normIDs <- c(normIDs, newNormIDs)
        print(paste('ncells for', p, 'greater than 400'))
      }
    }
  }
  
  filtered <- x[,normIDs]
  print(myTable(filtered))
  
  return(filtered)
}

DefaultAssay(s1p.combined) <- "RNA"

s1p.combined$patient_time <- paste(s1p.combined$patient, s1p.combined$time) %>%
  str_replace_all(" ", "_")
s1p.combined$patient_time_lineage <- paste(s1p.combined$patient, s1p.combined$time, s1p.combined$lineageClusterIdents)

s1p.temp <- subset(s1p.combined, lineageClusterIdents == 'CD4+ T cells' |
                     lineageClusterIdents == 'CD8+ T cells' |
                     lineageClusterIdents == 'B cells' |
                     lineageClusterIdents == 'NK cells' |
                     lineageClusterIdents == 'Monocytes' |
                     lineageClusterIdents == 'CD16+ Monocytes')
s1p.temp$lineageClusterIdents <- factor(s1p.temp$lineageClusterIdents, levels = unique(s1p.temp$lineageClusterIdents))

# Equalize the number of cells in each cluster for each sample down to 400 cells
s1p.combined.equalNum <- subsample(s1p.temp)

# Equalize the number of cells of PLACEBO patients in each cluster for each sample
s1p.combined.equalNum.placebo <- subsample(subset(s1p.temp, treatment == 'Placebo'))

# Equalize the number of cells of BAF312 patients in each cluster for each sample
s1p.combined.equalNum.BAF312 <- subsample(subset(s1p.temp, treatment == 'BAF312'))

save(s1p.combined.equalNum, s1p.combined.equalNum.placebo, s1p.combined.equalNum.BAF312, file = 's1p.combined_equalNum_down_to_400_cells.RData')




