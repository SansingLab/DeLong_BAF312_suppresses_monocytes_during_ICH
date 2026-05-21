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

mywd <- paste0("[path]/", myTime(), " scDEG by metadata equal num exported Tables")
dir.create(mywd)
setwd(mywd)

myList <- list()

col1 <- c('GSEA based on DEGs from single cell gene expression, populations downsampled',
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

colnames(description) <- c('GSEA comparisons by metadata', '')
myList[['Description']] <- description

formatResults <- function(x) {
  df <- as.data.frame(x)
  df$leadingEdge <- sapply(df$leadingEdge, toString)
  colnames(df) <- c('Pathway', 'p-value', 'padj', 'log2err', 'ES', 'NES', 'Size', 'Leading Edge', 'Enrichment')
  return(df)
}


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

addGSEA <- function(cellType) {
  tempList <- list()
  for (c in 1:length(comparisons)) {
    for (t in timepoints) {
      df <- GSEA_res_Metadata_Hallmark[[paste(comparisons[c], cellType, t)]]$Results
      df$leadingEdge <- sapply(df$leadingEdge, toString)
      
      colnames(df) <- c('Pathway', 'p-value', 'padj', 'log2err', 'ES', 'NES', 'Size', 'Leading Edge', 'Enrichment')
      
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

# Load GSEA results from scripts "11 - GSEA by metadata equal num down to ... .R"
load("[path]/GSEA_scDEGs_equalNum_down_to_1600_cells_Metadata_Hallmark.RData")
myList <- c(myList, addGSEA('Monocytes'))

load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Metadata_Hallmark.RData")
names(GSEA_res_Metadata_Hallmark) <- str_sub(names(GSEA_res_Metadata_Hallmark), start = 2) # remove leading space from these names
myList <- c(myList, addGSEA('CD16+ Monocytes'))
myList <- c(myList, addGSEA('NK cells'))

load("[path]/GSEA_scDEGs_equalNum_down_to_800_cells_Metadata_Hallmark.RData")
names(GSEA_res_Metadata_Hallmark) <- str_sub(names(GSEA_res_Metadata_Hallmark), start = 2) # remove leading space from these names
myList <- c(myList, addGSEA('CD4+ T cells'))

load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Metadata_Hallmark.RData")
names(GSEA_res_Metadata_Hallmark) <- str_sub(names(GSEA_res_Metadata_Hallmark), start = 2) # remove leading space from these names
myList <- c(myList, addGSEA('CD8+ T cells'))
myList <- c(myList, addGSEA('B cells'))

write_xlsx(myList, path = 'Supp Table 9 - GSEA by Metadata.xlsx', format_headers = FALSE)







