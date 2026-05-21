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

mywd <- paste0("[path]/", myTime(), " by metadata equalNum exported Tables")
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
          'B cells',
          '',
          'Category',
          'mRS Day 90',
          'NIHSS Day 1',
          'ICH Volume Day 1',
          'Perihematomal Edema Volume Day 1',
          'Relative Perihematomal Edema Volume (PHE Vol / ICH Vol) Day 1',
          'ICH score')

col2 <- c('',
          '',
          'Downsampled size',
          '1,600',
          '400',
          '400',
          '800',
          '400',
          '400',
          '',
          'Comparison',
          'Below 4 vs above or equal to 4',
          'Above or equal to 10 vs below 10',
          'Above or equal to 30mL vs below 30mL',
          'Above or equal to 30mL vs below 30mL',
          'Above or equal to 1 vs below 1',
          'Above or equal to 2 vs below 2'
          )

description <- matrix(nrow = length(col2), ncol = 2)
description[,1] <- col1
description[,2] <- col2

description <- data.frame(' ' = col1,
                          ' ' = col2)

colnames(description) <- c('DEG comparisons by Metadata', '')
myList[['Description']] <- description

timepoints <- c('Day 1', 'Day 3', 'Day 7')

comparisons <- c('MRS_Day90_split_at_4 MRS_Day90_below_4 - MRS_Day90_above_or_equal_to_4',
                 'NIHSS_Day1_split_at_10 NIHSS_Day1_above_or_equal_to_10 - NIHSS_Day1_below_10',
                 'Hematoma_Vol_Day1_split_at_30 Hematoma_Vol_Day1_above_or_equal_to_30 - Hematoma_Vol_Day1_below_30',
                 'PHE_Day1_split_at_30 PHE_Day1_above_or_equal_to_30 - PHE_Day1_below_30',
                 'Rel_PHE_Day1_split_at_1 Rel_PHE_Day1_above_or_equal_to_1 - Rel_PHE_Day1_below_1',
                 'ICH_score_split_at_2 ICH_score_above_or_equal_to_2 - ICH_score_below_2')

tabNames <- c('mRS_D90',
              'NIHSS_D1',
              'Hem_Vol_D1',
              'PHE_Day_1',
              'Rel_PHE_D1',
              'ICH_score')

addDEGs <- function(pathRoot, suffix, cellType) {
  tempList <- list()
  for (c in 1:length(comparisons)) {
    for (t in timepoints) {
            df <- read.delim(paste0(pathRoot, comparisons[c], suffix, cellType, " ", t, ".txt"))
      colnames(df)[1] <- "Gene"
      
      shortTime <- case_when(t == 'Day 1' ~ 'D1',
                             t == 'Day 3' ~ 'D3',
                             t == 'Day 7' ~ 'D7')
      shortcellType <- case_when(cellType == 'Monocytes' ~ 'Mono',
                                 cellType == 'CD16+ Monocytes' ~ 'CD16_Mono',
                                 cellType == 'NK cells' ~ 'NK_cell',
                                 cellType == 'CD4+ T cells' ~ 'CD4',
                                 cellType == 'CD8+ T cells' ~ 'CD8',
                                 cellType == 'B cells' ~ 'B_cell'
                                 )
      tempList[[paste(tabNames[c], shortcellType, shortTime)]] <- df
    }
  }
  return(tempList)
}

# Make a list of DEGs generated in "08 - DEGs metadata EqualNum down to 400/800/1600.R" scripts)
myList <- c(myList, addDEGs("[path]", " equalNum down to 1600 cells ", "Monocytes"))
myList <- c(myList, addDEGs("[path]", " equalNum down to 400 cells ", "CD16+ Monocytes"))
myList <- c(myList, addDEGs("[path]", " equalNum down to 400 cells ", "NK cells"))
myList <- c(myList, addDEGs("[path]", " equalNum down to 800 cells ", "CD4+ T cells"))
myList <- c(myList, addDEGs("[path]", " equalNum down to 400 cells ", "CD8+ T cells"))
myList <- c(myList, addDEGs("[path]", " equalNum down to 400 cells ", "B cells"))

write_xlsx(myList, path = 'Supp Table 8 - DEG by Metadata.xlsx', format_headers = FALSE)







