library(dplyr)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(stringr)
library(scales) # to control ggplot color palette
library(writexl)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]/", myTime(), " BAF312 by time equalNum exported Tables")
dir.create(mywd)
setwd(mywd)

myList <- list()

# Create columns for "Description" tab of exported table
col1 <- c('DEGs from single cell gene expression, populations downsampled',
          '',
          'Population',
          'Monocytes',
          'CD16+ Monocytes',
          'NK cells',
          'CD4+ T cells',
          'CD8+ T cells',
          'B cells')

col2 <- c('',
          '',
          'Downsampled size',
          '1,600',
          '400',
          '400',
          '800',
          '400',
          '400')

description <- matrix(nrow = 9, ncol = 2)
description[,1] <- col1
description[,2] <- col2

description <- data.frame(' ' = col1,
                          ' ' = col2)

colnames(description) <- c('DEGs of BAF312-treated patients, comparisons by Day', '')
myList[['Description']] <- description

timepoints <- c('Day 3 vs Day 1', 'Day 7 vs Day 3', 'Day 7 vs Day 1')

addDEGs <- function(pathRoot, cellType) {
  tempList <- list()
  for (t in timepoints) {
    df <- read.delim(paste0(pathRoot, cellType, " ", t, ".txt"))
    colnames(df)[1] <- "Gene"
    shortTime <- case_when(t == 'Day 3 vs Day 1' ~ 'D3vD1',
                           t == 'Day 7 vs Day 3' ~ 'D7vD3',
                           t == 'Day 7 vs Day 1' ~ 'D7vD1')
    tempList[[paste('BAF312', cellType, shortTime)]] <- df
  }
  return(tempList)
}

# Make a list of DEGs generated in "08 - DEGs BAF312-treated EqualNum down to 400/800/1600.R" scripts)
myList <- c(myList, addDEGs("[path]/DEGs BAF312-treated equalNum down to 1600 cells ", "Monocytes"))
myList <- c(myList, addDEGs("[path]/DEGs BAF312-treated equalNum down to 400 cells ", "CD16+ Monocytes"))
myList <- c(myList, addDEGs("[path]/DEGs BAF312-treated equalNum down to 400 cells ", "NK cells"))
myList <- c(myList, addDEGs("[path]/DEGs BAF312-treated equalNum down to 800 cells ", "CD4+ T cells"))
myList <- c(myList, addDEGs("[path]/DEGs BAF312-treated equalNum down to 400 cells ", "CD8+ T cells"))
myList <- c(myList, addDEGs("[path]/DEGs BAF312-treated equalNum down to 400 cells ", "B cells"))

write_xlsx(myList, path = 'Supp Table 4 - DEG BAF312 by Day.xlsx', format_headers = FALSE)







