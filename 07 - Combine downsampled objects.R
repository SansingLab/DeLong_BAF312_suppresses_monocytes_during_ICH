library(Seurat)
library(cowplot)
library(tidyverse)
library(RColorBrewer)
library(patchwork) # for wrapping featureplots
library(viridis)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]/", myTime())
dir.create(mywd)
setwd(mywd)

# Create new seurat object from downsampled objects
load("[path]/s1p.combined_subclustered.RData")

# Load downsampled data from Scripts 06 - Equalize cell number down to ... cells.R
load("[path]/s1p.combined_equalNum_down_to_1600_cells.RData")
colnames1600 <- colnames(subset(s1p.combined.equalNum, lineageClusterIdents == "Monocytes"))
load("[path]/s1p.combined_equalNum_down_to_800_cells.RData")
colnames800 <- colnames(subset(s1p.combined.equalNum, lineageClusterIdents == "CD4+ T cells"))
load("[path]/s1p.combined_equalNum_down_to_400_cells.RData")
colnames400 <- colnames(subset(s1p.combined.equalNum, lineageClusterIdents == "CD8+ T cells" | lineageClusterIdents == "CD16+ Monocytes" | 
                                 lineageClusterIdents == "NK cells" | lineageClusterIdents == "B cells"))

s1p.downsampled <- subset(s1p.combined, cells = c(colnames1600, colnames800, colnames400))


# Raw Age
s1p.downsampled$Age <- NA
s1p.downsampled$Age[s1p.downsampled$patient == "s1p2"] <- 46
s1p.downsampled$Age[s1p.downsampled$patient == "s1p3"] <- 48
s1p.downsampled$Age[s1p.downsampled$patient == "s1p4"] <- 73
s1p.downsampled$Age[s1p.downsampled$patient == "s1p5"] <- 70
s1p.downsampled$Age[s1p.downsampled$patient == "s1p6"] <- 69
s1p.downsampled$Age[s1p.downsampled$patient == "s1p7"] <- 81
s1p.downsampled$Age[s1p.downsampled$patient == "s1p8"] <- 64

# Binary Age
s1p.downsampled$Age_split_at_60 <- NA
s1p.downsampled$Age_split_at_60[s1p.downsampled$Age >= 60] <- "Age_above_or_equal_to_60"
s1p.downsampled$Age_split_at_60[s1p.downsampled$Age < 60] <- "Age_below_60"

# Raw BMI
s1p.downsampled$BMI <- NA
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p2"] <- 42.1
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p3"] <- 29.5 
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p4"] <- 30.5 
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p5"] <- 40.7 
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p6"] <- 34.2 
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p7"] <- 29.1 
s1p.downsampled$BMI[s1p.downsampled$patient == "s1p8"] <- 27.4 

# Binary BMI
s1p.downsampled$BMI_split_at_32 <- NA
s1p.downsampled$BMI_split_at_32[s1p.downsampled$BMI >= 32] <- "BMI_above_or_equal_to_32"
s1p.downsampled$BMI_split_at_32[s1p.downsampled$BMI < 32] <- "BMI_below_32"

# Race
s1p.downsampled$Race <- NA
s1p.downsampled$Race[s1p.downsampled$patient == "s1p2"] <- "Caucasian"
s1p.downsampled$Race[s1p.downsampled$patient == "s1p3"] <- "Black"
s1p.downsampled$Race[s1p.downsampled$patient == "s1p4"] <- "Caucasian"
s1p.downsampled$Race[s1p.downsampled$patient == "s1p5"] <- "Caucasian"
s1p.downsampled$Race[s1p.downsampled$patient == "s1p6"] <- "Caucasian"
s1p.downsampled$Race[s1p.downsampled$patient == "s1p7"] <- "Black"
s1p.downsampled$Race[s1p.downsampled$patient == "s1p8"] <- "Caucasian"

# Sex
s1p.downsampled$Sex <- NA
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p2"] <- "Female"
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p3"] <- "Female"
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p4"] <- "Female"
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p5"] <- "Female"
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p6"] <- "Female"
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p7"] <- "Female"
s1p.downsampled$Sex[s1p.downsampled$patient == "s1p8"] <- "Male"

# Raw MRS values
s1p.downsampled$MRS_Day1 <- NA
s1p.downsampled$MRS_Day15 <- NA
s1p.downsampled$MRS_Day30 <- NA
s1p.downsampled$MRS_Day90 <- NA

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p2"] <- 0 
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p2"] <- 3 
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p2"] <- 2
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p2"] <- 2

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p3"] <- 0
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p3"] <- 5
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p3"] <- 3 
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p3"] <- 2 

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p4"] <- 0 
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p4"] <- 4 
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p4"] <- 4 
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p4"] <- 4 

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p5"] <- 0
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p5"] <- 4
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p5"] <- 3
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p5"] <- 3

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p6"] <- 0
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p6"] <- 5 
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p6"] <- 4 
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p6"] <- 4 

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p7"] <- 1 
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p7"] <- 5 
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p7"] <- 4 
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p7"] <- 4 

s1p.downsampled$MRS_Day1[s1p.downsampled$patient == "s1p8"] <- 0 
s1p.downsampled$MRS_Day15[s1p.downsampled$patient == "s1p8"] <- 2 
s1p.downsampled$MRS_Day30[s1p.downsampled$patient == "s1p8"] <- 2 
s1p.downsampled$MRS_Day90[s1p.downsampled$patient == "s1p8"] <- 1  

# Binary MRS outcome
s1p.downsampled$MRS_Day30_split_at_3 <- NA
s1p.downsampled$MRS_Day90_split_at_3 <- NA

s1p.downsampled$MRS_Day30_split_at_3[s1p.downsampled$MRS_Day30 >= 3] <- "MRS_Day30_above_or_equal_to_3"
s1p.downsampled$MRS_Day30_split_at_3[s1p.downsampled$MRS_Day30 < 3] <- "MRS_Day30_below_3"
s1p.downsampled$MRS_Day90_split_at_3[s1p.downsampled$MRS_Day90 >= 3] <- "MRS_Day90_above_or_equal_to_3"
s1p.downsampled$MRS_Day90_split_at_3[s1p.downsampled$MRS_Day90 < 3] <- "MRS_Day90_below_3"

s1p.downsampled$MRS_Day15_split_at_4 <- NA
s1p.downsampled$MRS_Day30_split_at_4 <- NA
s1p.downsampled$MRS_Day90_split_at_4 <- NA

s1p.downsampled$MRS_Day15_split_at_4[s1p.downsampled$MRS_Day15 >= 4] <- "MRS_Day15_above_or_equal_to_4"
s1p.downsampled$MRS_Day15_split_at_4[s1p.downsampled$MRS_Day15 < 4] <- "MRS_Day15_below_4"
s1p.downsampled$MRS_Day30_split_at_4[s1p.downsampled$MRS_Day30 >= 4] <- "MRS_Day30_above_or_equal_to_4"
s1p.downsampled$MRS_Day30_split_at_4[s1p.downsampled$MRS_Day30 < 4] <- "MRS_Day30_below_4"
s1p.downsampled$MRS_Day90_split_at_4[s1p.downsampled$MRS_Day90 >= 4] <- "MRS_Day90_above_or_equal_to_4"
s1p.downsampled$MRS_Day90_split_at_4[s1p.downsampled$MRS_Day90 < 4] <- "MRS_Day90_below_4"

s1p.downsampled$MRS_Day15_split_at_5 <- NA

s1p.downsampled$MRS_Day15_split_at_5[s1p.downsampled$MRS_Day15 >= 5] <- "MRS_Day15_above_or_equal_to_5"
s1p.downsampled$MRS_Day15_split_at_5[s1p.downsampled$MRS_Day15 < 5] <- "MRS_Day15_below_5"


# MRS change values
s1p.downsampled$MRS_change_from_d15_Day30 <- NA
s1p.downsampled$MRS_change_from_d15_Day90 <- NA

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p2"] <- -1
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p2"] <- -1

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p3"] <- -2 
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p3"] <- -3 

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p4"] <- 0 
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p4"] <- 0 

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p5"] <- -1
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p5"] <- -1

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p6"] <- -1 
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p6"] <- -1 

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p7"] <- -1 
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p7"] <- -1 

s1p.downsampled$MRS_change_from_d15_Day30[s1p.downsampled$patient == "s1p8"] <- 0 
s1p.downsampled$MRS_change_from_d15_Day90[s1p.downsampled$patient == "s1p8"] <- -1  


# Raw NIHSS values
s1p.downsampled$NIHSS_Day1 <- NA
s1p.downsampled$NIHSS_Day90 <- NA

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p2"] <- 14
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p2"] <- 2

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p3"] <- 16
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p3"] <- 4 

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p4"] <- 19
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p4"] <- 9 

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p5"] <- 3
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p5"] <- 2

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p6"] <- 10
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p6"] <- NA 

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p7"] <- 8
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p7"] <- NA

s1p.downsampled$NIHSS_Day1[s1p.downsampled$patient == "s1p8"] <- 4
s1p.downsampled$NIHSS_Day90[s1p.downsampled$patient == "s1p8"] <- NA   


# Binary NIHSS outcome
s1p.downsampled$NIHSS_Day1_split_at_10 <- NA
s1p.downsampled$NIHSS_Day90_split_at_4 <- NA

s1p.downsampled$NIHSS_Day1_split_at_10[s1p.downsampled$NIHSS_Day1 >= 10] <- "NIHSS_Day1_above_or_equal_to_10"
s1p.downsampled$NIHSS_Day1_split_at_10[s1p.downsampled$NIHSS_Day1 < 10] <- "NIHSS_Day1_below_10"
s1p.downsampled$NIHSS_Day90_split_at_4[s1p.downsampled$NIHSS_Day90 >= 4] <- "NIHSS_Day90_above_or_equal_to_4"
s1p.downsampled$NIHSS_Day90_split_at_4[s1p.downsampled$NIHSS_Day90 < 4] <- "NIHSS_Day90_below_4"

# NIHSS change values
s1p.downsampled$NIHSS_change_from_d1_Day90 <- NA

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p2"] <- -12

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p3"] <- -12

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p4"] <- -10

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p5"] <- -1

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p6"] <- NA 

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p7"] <- NA 

s1p.downsampled$NIHSS_change_from_d1_Day90[s1p.downsampled$patient == "s1p8"] <- NA  


# Raw absolute perihematomal edema volume (mL)
s1p.downsampled$PHE_Day1 <- NA
s1p.downsampled$PHE_Day8 <- NA
s1p.downsampled$PHE_Day15 <- NA

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p2"] <- 8.94
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p2"] <- 16.95
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p2"] <- 6.07

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p3"] <- 20.83
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p3"] <- 44.5
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p3"] <- 54.05

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p4"] <- 29.94
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p4"] <- 70.17
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p4"] <- 94.33

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p5"] <- 39.16
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p5"] <- 46.94
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p5"] <- 38.44

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p6"] <- 12.93
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p6"] <- 35.48
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p6"] <- 61.96

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p7"] <- 138.68 
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p7"] <- 170.75
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p7"] <- 176.85

s1p.downsampled$PHE_Day1[s1p.downsampled$patient == "s1p8"] <- 38.66
s1p.downsampled$PHE_Day8[s1p.downsampled$patient == "s1p8"] <- 33.12
s1p.downsampled$PHE_Day15[s1p.downsampled$patient == "s1p8"] <- 38.61

# Binary absolute perihematomal edema volume
s1p.downsampled$PHE_Day1_split_at_30 <- NA
s1p.downsampled$PHE_Day8_split_at_40 <- NA
s1p.downsampled$PHE_Day15_split_at_40 <- NA

s1p.downsampled$PHE_Day1_split_at_30[s1p.downsampled$PHE_Day1 >= 30] <- "PHE_Day1_above_or_equal_to_30"
s1p.downsampled$PHE_Day1_split_at_30[s1p.downsampled$PHE_Day1 < 30] <- "PHE_Day1_below_30"
s1p.downsampled$PHE_Day8_split_at_40[s1p.downsampled$PHE_Day8 >= 40] <- "PHE_Day8_above_or_equal_to_40"
s1p.downsampled$PHE_Day8_split_at_40[s1p.downsampled$PHE_Day8 < 40] <- "PHE_Day8_below_40"
s1p.downsampled$PHE_Day15_split_at_40[s1p.downsampled$PHE_Day15 >= 40] <- "PHE_Day15_above_or_equal_to_40"
s1p.downsampled$PHE_Day15_split_at_40[s1p.downsampled$PHE_Day15 < 40] <- "PHE_Day15_below_40"


# Raw hematoma volume (mL)
s1p.downsampled$Hematoma_Vol_Day1 <- NA
s1p.downsampled$Hematoma_Vol_Day8 <- NA
s1p.downsampled$Hematoma_Vol_Day15 <- NA

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p2"] <- 7.28
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p2"] <- 9.48
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p2"] <- 1.9

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p3"] <- 21.36
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p3"] <- 21.9
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p3"] <- 13.63

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p4"] <- 45.9
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p4"] <- 42.04
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p4"] <- 26.5

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p5"] <- 34.32
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p5"] <- 4.28
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p5"] <- 2.52

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p6"] <- 24.48
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p6"] <- 20.07
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p6"] <- 9.38

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p7"] <- 48.47
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p7"] <- 54.01
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p7"] <- 30.28

s1p.downsampled$Hematoma_Vol_Day1[s1p.downsampled$patient == "s1p8"] <- 26.63
s1p.downsampled$Hematoma_Vol_Day8[s1p.downsampled$patient == "s1p8"] <- 16.67
s1p.downsampled$Hematoma_Vol_Day15[s1p.downsampled$patient == "s1p8"] <- 5.56

# Binary hematoma volume outcome
s1p.downsampled$Hematoma_Vol_Day1_split_at_30 <- NA
s1p.downsampled$Hematoma_Vol_Day8_split_at_30 <- NA
s1p.downsampled$Hematoma_Vol_Day15_split_at_20 <- NA

s1p.downsampled$Hematoma_Vol_Day1_split_at_30[s1p.downsampled$Hematoma_Vol_Day1 >= 30] <- "Hematoma_Vol_Day1_above_or_equal_to_30"
s1p.downsampled$Hematoma_Vol_Day1_split_at_30[s1p.downsampled$Hematoma_Vol_Day1 < 30] <- "Hematoma_Vol_Day1_below_30"
s1p.downsampled$Hematoma_Vol_Day8_split_at_30[s1p.downsampled$Hematoma_Vol_Day8 >= 30] <- "Hematoma_Vol_Day8_above_or_equal_to_30"
s1p.downsampled$Hematoma_Vol_Day8_split_at_30[s1p.downsampled$Hematoma_Vol_Day8 < 30] <- "Hematoma_Vol_Day8_below_30"
s1p.downsampled$Hematoma_Vol_Day15_split_at_20[s1p.downsampled$Hematoma_Vol_Day15 >= 20] <- "Hematoma_Vol_Day15_above_or_equal_to_20"
s1p.downsampled$Hematoma_Vol_Day15_split_at_20[s1p.downsampled$Hematoma_Vol_Day15 < 20] <- "Hematoma_Vol_Day15_below_20"


# Raw relative perihematomal edema volume (mL)
s1p.downsampled$Rel_PHE_Day1 <- NA
s1p.downsampled$Rel_PHE_Day8 <- NA
s1p.downsampled$Rel_PHE_Day15 <- NA

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p2"] <- 0.551
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p2"] <- 0.641
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p2"] <- 0.762

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p3"] <- 0.493
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p3"] <- 0.67
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p3"] <- 0.799

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p4"] <- 0.395
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p4"] <- 0.625
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p4"] <- 0.781

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p5"] <- 0.533
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p5"] <- 0.916
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p5"] <- 0.938

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p6"] <- 0.346
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p6"] <- 0.639
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p6"] <- 0.869

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p7"] <- 0.741
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p7"] <- 0.76
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p7"] <- 0.854

s1p.downsampled$Rel_PHE_Day1[s1p.downsampled$patient == "s1p8"] <- 0.592
s1p.downsampled$Rel_PHE_Day8[s1p.downsampled$patient == "s1p8"] <- 0.665
s1p.downsampled$Rel_PHE_Day15[s1p.downsampled$patient == "s1p8"] <- 0.874


# Binary relative perihematomal edema volume outcome
s1p.downsampled$Rel_PHE_Day1_split_at_0_5 <- NA
s1p.downsampled$Rel_PHE_Day8_split_at_0_65 <- NA
s1p.downsampled$Rel_PHE_Day15_split_at_0_8 <- NA

s1p.downsampled$Rel_PHE_Day1_split_at_0_5[s1p.downsampled$Rel_PHE_Day1 >= 0.5] <- "Rel_PHE_Day1_above_or_equal_to_0_5"
s1p.downsampled$Rel_PHE_Day1_split_at_0_5[s1p.downsampled$Rel_PHE_Day1 < 0.5] <- "Rel_PHE_Day1_below_0_5"
s1p.downsampled$Rel_PHE_Day8_split_at_0_65[s1p.downsampled$Rel_PHE_Day8 >= 0.65] <- "Rel_PHE_Day8_above_or_equal_to_0_65"
s1p.downsampled$Rel_PHE_Day8_split_at_0_65[s1p.downsampled$Rel_PHE_Day8 < 0.65] <- "Rel_PHE_Day8_below_0_65"
s1p.downsampled$Rel_PHE_Day15_split_at_0_8[s1p.downsampled$Rel_PHE_Day15 >= 0.8] <- "Rel_PHE_Day15_above_or_equal_to_0_8"
s1p.downsampled$Rel_PHE_Day15_split_at_0_8[s1p.downsampled$Rel_PHE_Day15 < 0.8] <- "Rel_PHE_Day15_below_0_8"


# Hemorrhage location
s1p.downsampled$Hematoma_Location <- NA

s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p2"] <- "Basal_ganglia"
s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p3"] <- "Basal_ganglia"
s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p4"] <- "Basal_ganglia"
s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p5"] <- "Occipital_lobe"
s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p6"] <- "Basal_ganglia"
s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p7"] <- "Occipital_lobe" 
s1p.downsampled$Hematoma_Location[s1p.downsampled$patient == "s1p8"] <- "Occipital_lobe"


# ICH Score
s1p.downsampled$ICH_score <- NA

s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p2'] <- 0
s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p3'] <- 0
s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p4'] <- 2
s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p5'] <- 1
s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p6'] <- 1
s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p7'] <- 2
s1p.downsampled$ICH_score[s1p.downsampled$patient == 's1p8'] <- 0

# Binary ICH score
s1p.downsampled$ICH_score_split_at_1 <- NA
s1p.downsampled$ICH_score_split_at_1[s1p.downsampled$ICH_score >= 1] <- "ICH_score_above_or_equal_to_1"
s1p.downsampled$ICH_score_split_at_1[s1p.downsampled$ICH_score < 1] <- "ICH_score_below_1"

s1p.downsampled$ICH_score_split_at_2 <- NA
s1p.downsampled$ICH_score_split_at_2[s1p.downsampled$ICH_score >= 2] <- "ICH_score_above_or_equal_to_2"
s1p.downsampled$ICH_score_split_at_2[s1p.downsampled$ICH_score < 2] <- "ICH_score_below_2"

save(s1p.downsampled, file = "s1p.downsampled.RData")
