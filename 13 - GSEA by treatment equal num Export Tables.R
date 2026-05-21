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

mywd <- paste0("[path]", myTime(), " scDEG by treatment equal num exported Tables")
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

colnames(description) <- c('GSEA comparisons by treatment (BAF312 vs Placebo)', '')
myList[['Description']] <- description

formatResults <- function(x) {
  df <- as.data.frame(x)
  df$leadingEdge <- sapply(df$leadingEdge, toString)
  colnames(df) <- c('Pathway', 'p-value', 'padj', 'log2err', 'ES', 'NES', 'Size', 'Leading Edge', 'Enrichment')
  return(df)
}

# Load GSEA results from scripts "11 - GSEA by treatment equal num down to ... .R"
load("[path]/GSEA_scDEGs_equalNum_down_to_1600_cells_Treatment_Hallmark.RData")
myList[['Monocytes Day 1']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment Monocytes Day 1`$Results)
myList[['Monocytes Day 3']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment Monocytes Day 3`$Results)
myList[['Monocytes Day 7']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment Monocytes Day 7`$Results)

load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Treatment_Hallmark.RData")
myList[['CD16+ Monocytes Day 1']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD16+ Monocytes Day 1`$Results)
myList[['CD16+ Monocytes Day 3']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD16+ Monocytes Day 3`$Results)
myList[['CD16+ Monocytes Day 7']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD16+ Monocytes Day 7`$Results)

myList[['NK cells Day 1']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment NK cells Day 1`$Results)
myList[['NK cells Day 3']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment NK cells Day 3`$Results)
myList[['NK cells Day 7']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment NK cells Day 7`$Results)

load("[path]/GSEA_scDEGs_equalNum_down_to_800_cells_Treatment_Hallmark.RData")
myList[['CD4+ T cells Day 1']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD4+ T cells Day 1`$Results)
myList[['CD4+ T cells Day 3']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD4+ T cells Day 3`$Results)
myList[['CD4+ T cells Day 7']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD4+ T cells Day 7`$Results)

load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Treatment_Hallmark.RData")
myList[['CD8+ T cells Day 1']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD8+ T cells Day 1`$Results)
myList[['CD8+ T cells Day 3']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD8+ T cells Day 3`$Results)
myList[['CD8+ T cells Day 7']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment CD8+ T cells Day 7`$Results)

myList[['B cells Day 1']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment B cells Day 1`$Results)
myList[['B cells Day 3']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment B cells Day 3`$Results)
myList[['B cells Day 7']] <- formatResults(GSEA_res_Treatment_Hallmark$`Treatment B cells Day 7`$Results)

write_xlsx(myList, path = 'Supp Table 7 - GSEA by Treatment.xlsx', format_headers = FALSE)







