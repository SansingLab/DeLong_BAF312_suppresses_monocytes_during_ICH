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

mywd <- paste0("[path]/", myTime(), " individual points")
dir.create(mywd)
setwd(mywd)

# Load clustered data with metadata added from Script 03 - Add Patient Metadata
load("[path]/s1p.combined_metadata.RData")

patientList = c("s1p2", "s1p3", "s1p4", "s1p5", "s1p6", "s1p7", "s1p8")

patientListPlacebo = c("s1p5", "s1p6")
patientListBAF312 = c("s1p3", "s1p4", "s1p7", "s1p8")
patientListBAF312_interrupted = c("s1p2")

dayList = c("Day 1", "Day 3", "Day 7")

DefaultAssay(s1p.combined) <- "RNA"

s1p.combined <- NormalizeData(s1p.combined, assay = "RNA")

# Functions for making error bars
se <- function(x) sd(x)/sqrt(length(x))
meanminusse <- function(x) return (mean(x) - se(x))
meanplusse <- function(x) return (mean(x) + se(x))


# # # Create table and graphs for broadClusters # # #
s1p.combined <- SetIdent(s1p.combined, value = 'broadClusterIdents')

# Print freq table (# of total cells) -- used for freq of singlets tables
df <- data.frame(Patient = character(), Day = character(), Cluster = character(), Freq = double())
for (d in dayList) {
  s1p.day <- subset(s1p.combined, subset = time == d)
  for (p in patientList) {
    new_df <- as.data.frame(sort(prop.table(table(Idents(subset(s1p.day, subset = patient == p))))))
    
    colnames(new_df) <- c("Cluster", "Freq")
    rownames(new_df) <- NULL
    new_df$Cluster <- as.character(new_df$Cluster)
    new_df$Freq <- as.numeric(new_df$Freq)
    new_df$Patient <- p
    new_df$Day <- d
    df <- dplyr::union(df, new_df)
  }
}
df$Treatment[df$Patient %in% patientListPlacebo] <- "Placebo"
df$Treatment[df$Patient %in% patientListBAF312] <- "BAF312"
df$Treatment[df$Patient %in% patientListBAF312_interrupted] <- "BAF312_interrupted"

gData <- df
gData <- filter(gData, Treatment != "BAF312_interrupted")

gData$Percent <- 100 * gData$Freq
gData$Treatment <- as.factor(gData$Treatment)
gData$Treatment <- factor(gData$Treatment, levels = c("Placebo", "BAF312"))
levels(gData$Treatment)
gData$Cluster <- factor(gData$Cluster, levels = rev(levels(s1p.combined$broadClusterIdents)))
levels(gData$Cluster)

# Graph broadclusters ratio vs Day 1 in BAF312-treated patients, T cells and B cells only
# Add column for percent of Day 1
df_spread <- spread(df, Day, Freq) %>%
  filter(Cluster %in% c('CD4 Naive', 'CD4 TCM', 'CD4 TEM',
                        'Treg',
                        'CD8 Naive', 'CD8 TCM', 'CD8 TEM',
                        'B Naive', 'B Int', 'B Mem')) %>%
  replace(is.na(.), 0) %>% # deal with populations not present in the full df
  mutate(Day3v1_Ratio = `Day 3` / `Day 1`) %>%
  mutate(Day7v1_Ratio = `Day 7` / `Day 1`)

df_ratio <- gather(df_spread, 'Ratio_Type', 'Ratio', 7:8)

gData <- df_ratio
gData$Treatment <- as.factor(gData$Treatment)
gData$Treatment <- factor(gData$Treatment, levels = c("Placebo", "BAF312"))
levels(gData$Treatment)
gData$Cluster <- factor(gData$Cluster, levels = c('CD4 Naive', 'CD4 TCM', 'CD4 TEM',
                                                  'Treg',
                                                  'CD8 Naive', 'CD8 TCM', 'CD8 TEM',
                                                  'B Naive', 'B Int', 'B Mem'))
levels(gData$Cluster)

gData <- subset(gData, Treatment == "BAF312")

colors = setNames(c(hue_pal()(3)[1], hue_pal()(3)[2]),
                  c("Day3v1_Ratio", "Day7v1_Ratio"))

ggplot(gData, aes(x = Cluster, y = Ratio, group = Ratio_Type, fill = Ratio_Type)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  stat_summary(fun = mean, position = position_dodge(width = 0.8), geom = "bar", width = 0.7, alpha = 0.5) +
  stat_summary(fun = mean, position = position_dodge(width = 0.8), geom = "bar", width = 0.7, fill = NA, aes(color = Ratio_Type)) +
  stat_summary(fun = mean, position = position_dodge(width = 0.8), fun.min = "meanminusse", fun.max = "meanplusse",
               geom = "errorbar", width = 0.3, alpha = 0.5, aes(color = Ratio_Type)) +
  geom_point(position = position_dodge(width = 0.8), shape = 16, aes(color = Ratio_Type)) +
  geom_hline(yintercept = 1, linetype = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))
ggsave2(filename = paste0("Ratio broadClusters T cells B cells BAF312-treated D3vD1 D7vD1 single points + error bars ", myTime(), ".pdf"), width = 6, height = 4)

write_delim(gData, paste0("Ratio broadClusters T cells B cells BAF312-treated D3vD1 D7vD1 single points + error bars table ", myTime (), ".txt"), delim = "\t")

gData_temp <- gData %>%
  group_by(Treatment, Cluster, Ratio_Type) %>%
  summarize(Average = mean(Ratio))
write_delim(gData_temp, paste0("Ratio broadClusters T cells B cells BAF312-treated D3vD1 D7vD1 single points + error bars table - average ratios ", myTime (), ".txt"), delim = "\t")


# # # Create table and graphs for lineageClusters # # #

s1p.combined <- SetIdent(s1p.combined, value = 'lineageClusterIdents')

# Print freq table (# of total cells) -- used for freq of singlets tables
df <- data.frame(Patient = character(), Day = character(), Cluster = character(), Freq = double())
for (d in dayList) {
  s1p.day <- subset(s1p.combined, subset = time == d)
  for (p in patientList) {
    new_df <- as.data.frame(sort(prop.table(table(Idents(subset(s1p.day, subset = patient == p))))))
    
    colnames(new_df) <- c("Cluster", "Freq")
    rownames(new_df) <- NULL
    new_df$Cluster <- as.character(new_df$Cluster)
    new_df$Freq <- as.numeric(new_df$Freq)
    new_df$Patient <- p
    new_df$Day <- d
    df <- dplyr::union(df, new_df)
  }
}

df$Treatment[df$Patient %in% patientListPlacebo] <- "Placebo"
df$Treatment[df$Patient %in% patientListBAF312] <- "BAF312"
df$Treatment[df$Patient %in% patientListBAF312_interrupted] <- "BAF312_interrupted"

df$Percent <- 100 * df$Freq

gData <- df
gData <- filter(gData, Treatment != "BAF312_interrupted")
gData$Treatment <- as.factor(gData$Treatment)
gData$Treatment <- factor(gData$Treatment, levels = c("Placebo", "BAF312"))
levels(gData$Treatment)
gData$Cluster <- factor(gData$Cluster, levels = rev(levels(s1p.combined$lineageClusterIdents)))
levels(gData$Cluster)


# Plot lymphocytes and myeloid cells in placebo-treated vs BAF312-treated patients

# Add column to indicate myeloid or lymphoid lineage
lymphoidList <- c('B cells', 'CD4+ T cells', 'CD8+ T cells', 'dnT cells', 'MAIT', 
                  'NK cells', 'pDC', 'Plasmablasts')
myeloidList <- c('CD16+ Monocytes', 'cDC1', 'cDC2', 'Monocytes')

adaptiveList <- c('B cells', 'CD4+ T cells', 'CD8+ T cells', 'dnT cells', 
                  'pDC', 'Plasmablasts')
innateList <- c('CD16+ Monocytes', 'cDC1', 'cDC2', 'Monocytes', 'NK cells', 'MAIT')

gData$Lineage <- case_when(gData$Cluster %in% lymphoidList ~ "Lymphoid",
                           gData$Cluster %in% myeloidList ~ "Myeloid",
                           TRUE ~ "Other")
gData$Innate <- case_when(gData$Cluster %in% innateList ~ "Innate",
                          gData$Cluster %in% adaptiveList ~ "Adaptive")


pal.lymphoid <- carto_pal(12, 'Bold')
pal.myeloid <- rev(carto_pal(12, 'Bold'))[-1]

# Facet by myeloid vs lymphoid
filter(gData, Cluster %in% c('Monocytes', 'CD16+ Monocytes', 'cDC1', 'cDC2', 'CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells')) %>%
  ggplot(aes(x = Day, y = Percent, group = Cluster, color = Cluster)) +
  geom_point(alpha = 0.5, shape = 16) +
  stat_summary(fun = mean, geom = "line", size = 1, alpha = 0.7, lineend = "round") +
  stat_summary(fun = mean, fun.min = "meanminusse", fun.max = "meanplusse",
               geom = "errorbar", width = 0.1, alpha = 0.4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text.x = element_blank()) +
  scale_color_manual(values = pal.lymphoid) +
  facet_wrap(Treatment ~ Innate, scales = "free", ncol = 2)
ggsave2(filename = paste0("Percent (major) lineageClusters by lineage and treatment individual points error bars ", myTime(), ".pdf"), width = 5, height = 4)

write_delim(gData, paste0("Percent (major) lineageClusters by lineage and treatment individual points error bars table ", myTime (), ".txt"), delim = "\t")

# Stacked barchart of lineage ratios split by patient
gData$Patient <- factor(gData$Patient, levels = c('s1p5', 's1p6', 's1p2', 's1p3', 's1p4', 's1p7', 's1p8'))
filter(gData, Cluster %in% c('Monocytes', 'CD16+ Monocytes', 'cDC1', 'cDC2', 'CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells')) %>%
  ggplot(aes(x = Patient, y = Percent, group = Cluster, fill = Cluster)) +
  geom_col(alpha = 0.7, color = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text.x = element_blank()) +
  scale_color_manual(values = pal.lymphoid) +
  scale_fill_manual(values = pal.lymphoid) +
  facet_wrap(~ Day, scales = "fixed", ncol = 3)
ggsave2(filename = paste0("Percent (major) lineageClusters by patient and day stacked bar charts ", myTime(), ".pdf"), width = 8, height = 3)


save.image(file = "Proportions.RData")
load('[path]/Proportions.RData')

