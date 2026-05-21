library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(sctransform)

# Palette libraries
library(RColorBrewer)
library(colorspace)
library(ggthemes)
library(scales)
library(rcartocolor)

library(ggh4x) # for splitting graph axes

# Function to return current date and time in a filename-friendly format
myTime <- function() {
        t <- gsub(":", "", Sys.time())
        t <- gsub(" ", "_", t)
        t <- gsub("\\.\\d+", "", t)
        return(t)
}

start_time <- Sys.time()

mywd <- paste0("[path]/", myTime())
dir.create(mywd)
setwd(mywd)

# Load clustered data with metadata added from Script 03 - Add Patient Metadata
load("[path]/s1p.combined_metadata.RData")

patientList = c("s1p2", "s1p3", "s1p4", "s1p5", "s1p6", "s1p7", "s1p8")
dayList = c("Day 1", "Day 3", "Day 7")

DefaultAssay(s1p.combined) <- "RNA"

s1p.combined <- NormalizeData(s1p.combined, assay = "RNA")

# Functions for making error bars
se <- function(x) sd(x)/sqrt(length(x))
meanminusse <- function(x) return (mean(x) - se(x))
meanplusse <- function(x) return (mean(x) + se(x))


# # # Create table and graphs for subclusters # # #
s1p.combined <- SetIdent(s1p.combined, value = 'lineageClusterIdents')

# Create freq table (# of total cells) -- used for freq of singlets tables
df <- data.frame(Patient = character(), Treatment = character(), Day = character(), Cluster = character(), Number = double())
for (d in dayList) {
  s1p.day <- subset(s1p.combined, subset = time == d)
  for (p in patientList) {
    s1p.day.patient <- subset(s1p.day, subset = patient == p)
    new_df <- as.data.frame(sort(table(Idents(s1p.day.patient))))
    
    colnames(new_df) <- c("Cluster", "Number")
    rownames(new_df) <- NULL
    new_df$Cluster <- as.character(new_df$Cluster)
    new_df$Number <- as.numeric(new_df$Number)
    new_df$Treatment <- unique(s1p.day.patient$treatment)
    new_df$Patient <- p
    new_df$Day <- d
    df <- dplyr::union(df, new_df)
  }
}

# Graph number of cells in each cluster
gData <- df

gData$Treatment <- as.factor(gData$Treatment)
gData$Treatment <- factor(gData$Treatment, levels = c("Placebo", "BAF312"))
levels(gData$Treatment)
gData$Cluster <- factor(gData$Cluster, levels = rev(levels(s1p.combined$lineageClusterIdents)))
levels(gData$Cluster)
gData$Patient <- str_replace(gData$Patient, "s1p", "Patient ")
gData <- filter(gData, Cluster %in% c('Monocytes', 'CD16+ Monocytes', 'NK cells', 'CD4+ T cells', 'CD8 T cells', 'B cells'))
gData$hline <- case_when(gData$Cluster == 'Monocytes' ~ 1600,
                         gData$Cluster == 'CD4+ T cells' ~ 800,
                         .default = 400)


ggplot(gData, aes(x = Day, y = Number, fill = Patient)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
  geom_hline(aes(yintercept = hline), size = 0.5, color = 'gray60', linetype = 'dashed') +
  facet_grid(rows = 'Cluster', scales = 'free')
ggsave(filename = paste0('Total cell number by cluster and day ', myTime(), '.pdf'), height = 7, width = 7)

write_delim(spread(gData, key = Day, value = Number), file = paste0('Total cell number by cluster and day table ', myTime(), '.txt'), delim = '\t')

