library(Seurat)
library(cowplot)
library(tidyverse)
library(RColorBrewer)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]", myTime(), " by treatment")
dir.create(mywd)
setwd(mywd)

timepoints <- c("Day 1", "Day 3", "Day 7")

# Function to load dataset
loadData <- function(cluster) {
  cluster <- str_replace(cluster, "\\+", "\\\\\\+")
  df <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(df) <- c('Gene', 
                    'avg_log2FC',
                    'p_val_adj',
                    'timepoint')
  
  df <- data.frame(Gene = character(), 
                   avg_log2FC = numeric(),
                   p_val_adj = numeric(),
                   timepoint = character())
  
  file_list <- list.files(path = DEGdir)
  file_list <- file_list[grep("txt$", file_list)]
  file_list <- file_list[grep("^DEGs by treatment equalNum down to", file_list)] # Filter on comparisons by treatment
  file_list_lineage <- file_list[grep(cluster, file_list)]
  
  if (length(file_list_lineage) > 0) { # Make sure there is at least one file for this lineage
    if (cluster == "Monocytes")
      file_list_lineage <- file_list_lineage[grep("CD16", file_list_lineage, invert = TRUE)]
    
    for (t in timepoints) {
      temp_data <- data.frame(matrix(data = NA, nrow = 0, ncol = 0))
      if (length(file_list_lineage[grep(t, file_list_lineage)]) > 0) { # Make sure there's a file for this timepoint
        temp_data <- data.frame(matrix(data = NA, nrow = 0, ncol = 0))
        temp_data <- read.delim(paste0(DEGdir, "/", file_list_lineage[grep(t, file_list_lineage)]), 
                                header = TRUE, sep = "\t", check.names = TRUE)
        colnames(temp_data)[1] <- 'Gene'
        
        temp_data <- temp_data[c(1, 3, 6)] # Retain only gene, avg_log2FC, p_val_adj
        temp_data$timepoint <- t
        df <- bind_rows(df, temp_data)
      }
    }
    
    # Fill in non-significant avg_log2FC with zeroes
    df_log2FC <- select(df, c(Gene, timepoint, avg_log2FC))
    df_log2FC <- spread(df_log2FC, key = timepoint, value = avg_log2FC)
    df_log2FC <- gather(df_log2FC, key = 'timepoint', value = 'avg_log2FC', 2:ncol(df_log2FC))
    df <- left_join(df_log2FC, df, by = c('Gene', 'timepoint'), suffix = c("", ".y")) %>%
      select(-ends_with(".y")) # join dfs and remove duplicate column
    df[is.na(df)] <- 0     # Convert NaN z-scores to zero
    
    # Order by log2FC
    df <- df[order(df$avg_log2FC),]
    
    # Add -log(adj p-value)
    df$log_p_val_adj <- ifelse(df$p_val_adj == 0, 10^-300, df$p_val_adj)
    df$log_p_val_adj <- -log10(df$log_p_val_adj)
    
    return(df)
  }
}

DotPlotSelectedGenesTreatment <- function(cluster, myFeatures, geneset = "select genes") {
  myData <- filter(df, str_detect(Gene, paste(paste0('^', myFeatures, '$'), collapse = "|")))
  
  myData$Gene <- factor(myData$Gene, levels = myFeatures)
  
  # Find width of longest canonical pathway label
  geneLabels <- as.character(myData$Gene)
  longestGeneLabel <- max(nchar(geneLabels))
  
  # Graphing code from https://www.biostars.org/p/359307/
  ggplot(myData, aes(x = timepoint, y = Gene, size = log_p_val_adj, color = avg_log2FC)) + 
    geom_point(alpha = 1) + 
    theme_classic() + 
    scale_color_steps2(low = "blue", mid = "white",  high = "red", midpoint = 0, space = "Lab") + 
    scale_size(range = c(2, 8)) +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
          axis.title = element_text(size = 15)) +
    labs(y = 'Gene', x = 'Timepoint',
         size = '-log(adj p-value)', color = 'avg_log2FC',
         title = paste(cluster, "-", geneset, "\n\n\n"))
  ggsave(paste0("Dotplot DEGs by treatment ", cluster, " scDEGs ", geneset, " ", myTime(),  ".pdf"), 
         height = (1.9 + .25 * length(myFeatures)), width = (4.5 + 0.055 * longestGeneLabel), device = cairo_pdf)
}

# # # # # # # # # # # # # # # # # # # # # #
# Downsampled to 1,600 cells (monocytes)  #
# # # # # # # # # # # # # # # # # # # # # #

# Load DEGs from script 08 - DEGs by treatment equal num down to 1600.R
DEGdir <- "[path]/2024-06-10_123151 by treatment equalNum down to 1600 cells"

# Monocytes
df <- loadData('Monocytes')
cluster <- 'Monocytes'

myFeatures <- rev(c(
  'TNFAIP3',
  'TNFAIP6',
  'NFKB1',
  'NFKBIA',
  'FOSB',
  'JUN',
  'CD83',
  'EGR1',
  'CXCL8',
  'CXCR4',
  'GBP4',
  'RIPK2',
  'RGS1',
  'CLU',
  'ABCA1'
))
  
DotPlotSelectedGenesTreatment('Monocytes', myFeatures)



