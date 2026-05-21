library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(dplyr)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  return(t)
}

# Set wd to the path filtered IPA pathways from script "17 - Filter IPA Upstream regulators"
basedir <- paste0("[path]/For IPA/UR by Treatment")
mywd <- paste0(basedir, "/Graphs/Cytokines ", myTime())
dir.create(mywd)
setwd(mywd)

# # # Run script to filter out qiagen copywrite before running this script # # #

subdir <- "/Filtered p<0_05 >4 molecules Cytokines"

timepoints <- c("Day 1", "Day 3", "Day 7")

# Function to load dataset
loadData <- function(cellType) {
  cellType <- str_replace(cellType, "\\+", "\\\\\\+")

  df <- data.frame(Upstream.Regulator = character(), 
                   B.H.corrected.p.value = numeric(), 
                   Bias.corrected.z.score = numeric(),
                   timepoint = character(),
                   Molecule_count = numeric())
  
  file_list <- list.files(path = paste0(basedir, subdir))
  file_list <- file_list[grep("^UR by treatment", file_list)] # Filter on comparisons by time
  
  file_list <- file_list[grep("txt$", file_list)]
  file_list_lineage <- file_list[grep(cellType, file_list)]
  
  if (cellType == "Monocytes")
    file_list_lineage <- file_list_lineage[grep("CD16", file_list_lineage, invert = TRUE)]
  
  for (t in timepoints) {
    if (length(file_list_lineage[grep(t, file_list_lineage)]) > 0) {
      temp_data <- read.delim(paste0(basedir, subdir, "/", file_list_lineage[grep(t, file_list_lineage)]), 
                              header = TRUE, sep = "\t", check.names = TRUE)
    
      temp_data <- temp_data[c(2, 6, 11, 14)] # Retain only p-value, z-score, and molecule count
      temp_data$timepoint <- t
      colnames(temp_data)[4] <- "Molecule_count"
      if (nrow(temp_data) > 0)
        df <- bind_rows(df, temp_data)
    }
  }
  
  # Fill in non-significant z-scores with zeroes
  df_zscores <- select(df, c(Upstream.Regulator, timepoint, Bias.corrected.z.score))
  df_zscores <- spread(df_zscores, key = timepoint, value = Bias.corrected.z.score)
  df_zscores <- gather(df_zscores, key = 'timepoint', value = 'Bias.corrected.z.score', 2:ncol(df_zscores))
  df <- left_join(df_zscores, df, by = c('Upstream.Regulator', 'timepoint'), suffix = c("", ".y")) %>%
    select(-ends_with(".y")) # join dfs and remove duplicate column
  df[is.na(df)] <- 0     # Convert NaN z-scores to zero
  
  # Order by z-score
  df <- df[order(df$Bias.corrected.z.score),]
  
  df$X.log.B.H.p.value. <- -log(df$B.H.corrected.p.value)
  
  return(df)
}

# Function to graph horizontal bar chart
graphBarChartHorizontal <- function(myData, cellType) {
  if (nrow(myData) > 0) { # Make sure there's data in the dataframe
    
    myData <- filter(myData, str_detect(Upstream.Regulator, paste(paste0('^', myPathways_high_z, '$'), collapse = "|")))
    
    myData$Upstream.Regulator <- factor(myData$Upstream.Regulator, levels = rev(myPathways_high_z))
    
    myData$timepoint <- factor(myData$timepoint, levels = c('Day 1', 'Day 3', 'Day 7'))
    
    # Find width of longest canonical pathway label
    cpLabels <- as.character(myData$Upstream.Regulator)
    longestCPLabel <- max(nchar(cpLabels))
    
    ggplot(myData, aes(x = Upstream.Regulator, y = Bias.corrected.z.score, fill = timepoint, color = timepoint)) + 
      geom_bar(aes(alpha = X.log.B.H.p.value.), stat = "identity", linewidth = 0.4, width = 0.8, position = position_dodge(width = 0.9)) +
      scale_fill_discrete(breaks=c('Day 1', 'Day 3', 'Day 7')) +
      scale_color_discrete(breaks=c('Day 1', 'Day 3', 'Day 7')) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = .5),
            axis.title = element_text(size = 15)) +
      geom_text(aes(label = ifelse(Bias.corrected.z.score == 0, "", Molecule_count), group = timepoint), 
                color = "Black", angle = 90, fontface = "italic", size = 3, hjust = ifelse(myData$Bias.corrected.z.score >= 0, -.2, 1.3),
                position = position_dodge(width = 0.9)) +
      labs(y = 'Bias-corrected z-score', x = 'IPA Upstream Regulator',
           alpha = '-log(BH-adj p-value)', color = 'Timepoint', fill = 'Timepoint',
           title = paste(cellType, "- highest z-scores\n")) +
      scale_y_continuous(limits = c(ifelse(min(myData$Bias.corrected.z.score) < 0, min(myData$Bias.corrected.z.score)*1.5, 0), ifelse(max(myData$Bias.corrected.z.score) > 0, max(myData$Bias.corrected.z.score)*1.5, 0)))
    ggsave(paste0("Barchart Horizontal IPA UR by treatment highest z-scores ", cellType, " " , myTime(),  ".pdf"), 
           width = (3.3 + .35 * length(myPathways_high_z)), height = (3.5 + 0.055 * longestCPLabel), device = cairo_pdf) 
  }
}

# # # Monocytes # # #
df <- loadData("Monocytes")

# Filter dataset by z-score

df_high_z <- df %>%
  arrange(desc(abs(Bias.corrected.z.score))) %>%
  distinct(Upstream.Regulator, .keep_all = TRUE)

myPathways_high_z <- rev(df_high_z[1:10, 1])

graphBarChartHorizontal(df, "Monocytes")



