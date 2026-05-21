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

mywd <- paste0("[path]", myTime(), " scDEG by time BAF312 equal num exported Tables")
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

colnames(description) <- c('GSEA of BAF312-treated patients, comparisons by Day', '')
myList[['Description']] <- description

formatResults <- function(x) {
  df <- as.data.frame(x)
  df$leadingEdge <- sapply(df$leadingEdge, toString)
  colnames(df) <- c('Pathway', 'p-value', 'padj', 'log2err', 'ES', 'NES', 'Size', 'Leading Edge', 'Enrichment')
  return(df)
}

# Load GSEA results from scripts "11 - GSEA by time BAF312 equal num down to ... .R"
load("[path]/GSEA_scDEGs_equalNum_down_to_1600_cells_Time_BAF312_Hallmark.RData")
myList[['BAF312 Monocytes D3vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time Monocytes D3vD1`$Results)
myList[['BAF312 Monocytes D7vD3']] <- formatResults(GSEA_res_Time_Hallmark$`time Monocytes D7vD3`$Results)
myList[['BAF312 Monocytes D7vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time Monocytes D7vD1`$Results)

load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Time_BAF312_Hallmark.RData")
myList[['BAF312 CD16+ Monocytes D3vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time CD16+ Monocytes D3vD1`$Results)
myList[['BAF312 CD16+ Monocytes D7vD3']] <- formatResults(GSEA_res_Time_Hallmark$`time CD16+ Monocytes D7vD3`$Results)
myList[['BAF312 CD16+ Monocytes D7vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time CD16+ Monocytes D7vD1`$Results)

myList[['BAF312 NK cells D3vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time NK cells D3vD1`$Results)
myList[['BAF312 NK cells D7vD3']] <- formatResults(GSEA_res_Time_Hallmark$`time NK cells D7vD3`$Results)
myList[['BAF312 NK cells D7vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time NK cells D7vD1`$Results)

load("[path]/GSEA_scDEGs_equalNum_down_to_800_cells_Time_BAF312_Hallmark.RData")
myList[['BAF312 CD4+ T cells D3vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time CD4+ T cells D3vD1`$Results)
myList[['BAF312 CD4+ T cells D7vD3']] <- formatResults(GSEA_res_Time_Hallmark$`time CD4+ T cells D7vD3`$Results)
myList[['BAF312 CD4+ T cells D7vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time CD4+ T cells D7vD1`$Results)

load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Time_BAF312_Hallmark.RData")
myList[['BAF312 CD8+ T cells D3vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time CD8+ T cells D3vD1`$Results)
myList[['BAF312 CD8+ T cells D7vD3']] <- formatResults(GSEA_res_Time_Hallmark$`time CD8+ T cells D7vD3`$Results)
myList[['BAF312 CD8+ T cells D7vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time CD8+ T cells D7vD1`$Results)

myList[['BAF312 B cells D3vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time B cells D3vD1`$Results)
myList[['BAF312 B cells D7vD3']] <- formatResults(GSEA_res_Time_Hallmark$`time B cells D7vD3`$Results)
myList[['BAF312 B cells D7vD1']] <- formatResults(GSEA_res_Time_Hallmark$`time B cells D7vD1`$Results)

write_xlsx(myList, path = 'Supp Table 5 - GSEA BAF312 by Day.xlsx', format_headers = FALSE)







