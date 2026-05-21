library(dplyr)
library(tidyverse)
library(ggrepel)
library(stringr)
library(scales) # to control ggplot color palette

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]/", myTime(), " scDEG by time Placebo equal num summary by pathway")
dir.create(mywd)
setwd(mywd)

# Load GSEA results from scripts "11 - GSEA by time Placebo equal num down to ... .R"
load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Time_Placebo_Hallmark.RData")
GSEA_res_Time_Hallmark_400 <- GSEA_res_Time_Hallmark
load("[path]/GSEA_scDEGs_equalNum_down_to_800_cells_Time_Placebo_Hallmark.RData")
GSEA_res_Time_Hallmark_800 <- GSEA_res_Time_Hallmark
load("[path]/GSEA_scDEGs_equalNum_down_to_1600_cells_Time_Placebo_Hallmark_seed.RData")
GSEA_res_Time_Hallmark_1600 <- GSEA_res_Time_Hallmark
rm(GSEA_res_Time_Hallmark)

# Create new list with appropriate downsampling of each population
GSEA_res_Time_combined <- GSEA_res_Time_Hallmark_400 # Start with 400, since most samples are downsampled to 400

CD4s <- c('time CD4+ T cells D3vD1', 'time CD4+ T cells D7vD1', 'time CD4+ T cells D7vD3')
for (pop in CD4s) {
  GSEA_res_Time_combined[[pop]] <- GSEA_res_Time_Hallmark_800[[pop]]
}

Monocytes <- c('time Monocytes D3vD1', 'time Monocytes D7vD1', 'time Monocytes D7vD3')
for (pop in Monocytes) {
  GSEA_res_Time_combined[[pop]] <- GSEA_res_Time_Hallmark_1600[[pop]]
}


plotSummary <- function(p) {
  df <- data.frame()
  #myPathway <- 'INTERFERON GAMMA RESPONSE'
  myPathway <- p
  for (c in names(GSEA_res_Time_combined)) {
    newrow <- filter(GSEA_res_Time_combined[[c]]$Results, pathway == myPathway)
    if (nrow(newrow)  == 0) { # if this pathway isn't enriched in this population
      newrow <- data.frame(pathway = myPathway,
                           pval = 1,
                           padj = 1,
                           log2err = 0,
                           ES = 0,
                           NES = 0,
                           size = 0,
                           leadingEdge = "blank",
                           Enrichment = "None")
    }
    newrow$comparison <- c
    newrow$population <- str_sub(c, start = 6, end = -7)
    newrow$time <- str_sub(c, start = -5)
    newrow$log_p <- -log10(newrow$padj)
    newrow <- select(newrow, -leadingEdge) # this list object causes problems when binding empty rows
    if (nrow(df) == 0) {
      df <- newrow
    } else {
      df <- bind_rows(df, newrow)
    }
  }
  
  df$population <- factor(df$population, levels = rev(c('Monocytes', 'CD16+ Monocytes', 'CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells')))
  
  colors = setNames(c("red", "white", "blue"),
                    c("Up-regulated", "None", "Down-regulated"))
   
  ggplot(filter(df, time != 'D7vD3'), aes(time, population)) +
    theme_classic() +
    geom_point(aes(color = NES, alpha = log_p), shape = 15, size = 6) +
    geom_point(shape = 0, size = 6) +
    scale_color_gradient2(midpoint=0, low="blue", mid="white",
                          high="red", space ="Lab" ) +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
          axis.title = element_text(size = 15),
          title = element_text(size = 10)) +
    labs(y = paste('GSEA', myPathway), x = 'NES',
         alpha = '-log(p-val)', color = 'NES',
         title = paste0("GSEA Placebo\nby Time\nscDEGs equalNum summary\n", myPathway, "\n\n\n\n\n"))
  ggsave(paste0("Dotplot GSEA scDEGs NoNPerm Placebo by Time equalNum summary min_pct=0_1 ", myPathway, " ", myTime(),  ".pdf"), 
         height = 4, width = 3.5, device = cairo_pdf)
}

pathwayList <- c('TNFA SIGNALING VIA NFKB', 'INTERFERON ALPHA RESPONSE', 'INTERFERON GAMMA RESPONSE', 
                 'CHOLESTEROL HOMEOSTASIS', 'IL2 STAT5 SIGNALING', 'INFLAMMATORY RESPONSE', 'HYPOXIA')

for (path in pathwayList){
  plotSummary(path)
}

