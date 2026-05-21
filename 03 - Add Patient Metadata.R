library(Seurat)

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

# Load clustered data from Script 02 - Cluster.R
load("[path]/s1p.combined_subclustered.RData")

# Raw Age
s1p.combined$Age <- NA
s1p.combined$Age[s1p.combined$patient == "s1p2"] <- 46
s1p.combined$Age[s1p.combined$patient == "s1p3"] <- 48
s1p.combined$Age[s1p.combined$patient == "s1p4"] <- 73
s1p.combined$Age[s1p.combined$patient == "s1p5"] <- 70
s1p.combined$Age[s1p.combined$patient == "s1p6"] <- 69
s1p.combined$Age[s1p.combined$patient == "s1p7"] <- 81
s1p.combined$Age[s1p.combined$patient == "s1p8"] <- 64

# Binary Age
s1p.combined$Age_split_at_60 <- NA
s1p.combined$Age_split_at_60[s1p.combined$Age >= 60] <- "Age_above_or_equal_to_60"
s1p.combined$Age_split_at_60[s1p.combined$Age < 60] <- "Age_below_60"

# Raw BMI
s1p.combined$BMI <- NA
s1p.combined$BMI[s1p.combined$patient == "s1p2"] <- 42.1
s1p.combined$BMI[s1p.combined$patient == "s1p3"] <- 29.5 
s1p.combined$BMI[s1p.combined$patient == "s1p4"] <- 30.5 
s1p.combined$BMI[s1p.combined$patient == "s1p5"] <- 40.7 
s1p.combined$BMI[s1p.combined$patient == "s1p6"] <- 34.2 
s1p.combined$BMI[s1p.combined$patient == "s1p7"] <- 29.1 
s1p.combined$BMI[s1p.combined$patient == "s1p8"] <- 27.4 

# Binary BMI
s1p.combined$BMI_split_at_32 <- NA
s1p.combined$BMI_split_at_32[s1p.combined$BMI >= 32] <- "BMI_above_or_equal_to_32"
s1p.combined$BMI_split_at_32[s1p.combined$BMI < 32] <- "BMI_below_32"

# Race
s1p.combined$Race <- NA
s1p.combined$Race[s1p.combined$patient == "s1p2"] <- "Caucasian"
s1p.combined$Race[s1p.combined$patient == "s1p3"] <- "Black"
s1p.combined$Race[s1p.combined$patient == "s1p4"] <- "Caucasian"
s1p.combined$Race[s1p.combined$patient == "s1p5"] <- "Caucasian"
s1p.combined$Race[s1p.combined$patient == "s1p6"] <- "Caucasian"
s1p.combined$Race[s1p.combined$patient == "s1p7"] <- "Black"
s1p.combined$Race[s1p.combined$patient == "s1p8"] <- "Caucasian"

# Sex
s1p.combined$Sex <- NA
s1p.combined$Sex[s1p.combined$patient == "s1p2"] <- "Female"
s1p.combined$Sex[s1p.combined$patient == "s1p3"] <- "Female"
s1p.combined$Sex[s1p.combined$patient == "s1p4"] <- "Female"
s1p.combined$Sex[s1p.combined$patient == "s1p5"] <- "Female"
s1p.combined$Sex[s1p.combined$patient == "s1p6"] <- "Female"
s1p.combined$Sex[s1p.combined$patient == "s1p7"] <- "Female"
s1p.combined$Sex[s1p.combined$patient == "s1p8"] <- "Male"

# Raw MRS values
s1p.combined$MRS_Day1 <- NA
s1p.combined$MRS_Day15 <- NA
s1p.combined$MRS_Day30 <- NA
s1p.combined$MRS_Day90 <- NA

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p2"] <- 0 
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p2"] <- 3 
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p2"] <- 2
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p2"] <- 2

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p3"] <- 0
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p3"] <- 5
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p3"] <- 3 
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p3"] <- 2 

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p4"] <- 0 
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p4"] <- 4 
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p4"] <- 4 
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p4"] <- 4 

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p5"] <- 0
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p5"] <- 4
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p5"] <- 3
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p5"] <- 3

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p6"] <- 0
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p6"] <- 5 
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p6"] <- 4 
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p6"] <- 4 

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p7"] <- 1 
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p7"] <- 5 
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p7"] <- 4 
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p7"] <- 4 

s1p.combined$MRS_Day1[s1p.combined$patient == "s1p8"] <- 0 
s1p.combined$MRS_Day15[s1p.combined$patient == "s1p8"] <- 2 
s1p.combined$MRS_Day30[s1p.combined$patient == "s1p8"] <- 2 
s1p.combined$MRS_Day90[s1p.combined$patient == "s1p8"] <- 1  

# Binary MRS outcome
s1p.combined$MRS_Day30_split_at_3 <- NA
s1p.combined$MRS_Day90_split_at_3 <- NA

s1p.combined$MRS_Day30_split_at_3[s1p.combined$MRS_Day30 >= 3] <- "MRS_Day30_above_or_equal_to_3"
s1p.combined$MRS_Day30_split_at_3[s1p.combined$MRS_Day30 < 3] <- "MRS_Day30_below_3"
s1p.combined$MRS_Day90_split_at_3[s1p.combined$MRS_Day90 >= 3] <- "MRS_Day90_above_or_equal_to_3"
s1p.combined$MRS_Day90_split_at_3[s1p.combined$MRS_Day90 < 3] <- "MRS_Day90_below_3"

s1p.combined$MRS_Day15_split_at_4 <- NA
s1p.combined$MRS_Day30_split_at_4 <- NA
s1p.combined$MRS_Day90_split_at_4 <- NA

s1p.combined$MRS_Day15_split_at_4[s1p.combined$MRS_Day15 >= 4] <- "MRS_Day15_above_or_equal_to_4"
s1p.combined$MRS_Day15_split_at_4[s1p.combined$MRS_Day15 < 4] <- "MRS_Day15_below_4"
s1p.combined$MRS_Day30_split_at_4[s1p.combined$MRS_Day30 >= 4] <- "MRS_Day30_above_or_equal_to_4"
s1p.combined$MRS_Day30_split_at_4[s1p.combined$MRS_Day30 < 4] <- "MRS_Day30_below_4"
s1p.combined$MRS_Day90_split_at_4[s1p.combined$MRS_Day90 >= 4] <- "MRS_Day90_above_or_equal_to_4"
s1p.combined$MRS_Day90_split_at_4[s1p.combined$MRS_Day90 < 4] <- "MRS_Day90_below_4"

s1p.combined$MRS_Day15_split_at_5 <- NA

s1p.combined$MRS_Day15_split_at_5[s1p.combined$MRS_Day15 >= 5] <- "MRS_Day15_above_or_equal_to_5"
s1p.combined$MRS_Day15_split_at_5[s1p.combined$MRS_Day15 < 5] <- "MRS_Day15_below_5"


# MRS change values
s1p.combined$MRS_change_from_d15_Day30 <- NA
s1p.combined$MRS_change_from_d15_Day90 <- NA

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p2"] <- -1
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p2"] <- -1

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p3"] <- -2 
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p3"] <- -3 

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p4"] <- 0 
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p4"] <- 0 

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p5"] <- -1
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p5"] <- -1

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p6"] <- -1 
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p6"] <- -1 

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p7"] <- -1 
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p7"] <- -1 

s1p.combined$MRS_change_from_d15_Day30[s1p.combined$patient == "s1p8"] <- 0 
s1p.combined$MRS_change_from_d15_Day90[s1p.combined$patient == "s1p8"] <- -1  


# Raw NIHSS values
s1p.combined$NIHSS_Day1 <- NA
s1p.combined$NIHSS_Day90 <- NA

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p2"] <- 14
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p2"] <- 2

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p3"] <- 16
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p3"] <- 4 

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p4"] <- 19
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p4"] <- 9 

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p5"] <- 3
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p5"] <- 2

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p6"] <- 10
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p6"] <- NA 

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p7"] <- 8
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p7"] <- NA

s1p.combined$NIHSS_Day1[s1p.combined$patient == "s1p8"] <- 4
s1p.combined$NIHSS_Day90[s1p.combined$patient == "s1p8"] <- NA   


# Binary NIHSS outcome
s1p.combined$NIHSS_Day1_split_at_10 <- NA
s1p.combined$NIHSS_Day90_split_at_4 <- NA

s1p.combined$NIHSS_Day1_split_at_10[s1p.combined$NIHSS_Day1 >= 10] <- "NIHSS_Day1_above_or_equal_to_10"
s1p.combined$NIHSS_Day1_split_at_10[s1p.combined$NIHSS_Day1 < 10] <- "NIHSS_Day1_below_10"
s1p.combined$NIHSS_Day90_split_at_4[s1p.combined$NIHSS_Day90 >= 4] <- "NIHSS_Day90_above_or_equal_to_4"
s1p.combined$NIHSS_Day90_split_at_4[s1p.combined$NIHSS_Day90 < 4] <- "NIHSS_Day90_below_4"

# NIHSS change values
s1p.combined$NIHSS_change_from_d1_Day90 <- NA

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p2"] <- -12

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p3"] <- -12

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p4"] <- -10

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p5"] <- -1

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p6"] <- NA 

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p7"] <- NA 

s1p.combined$NIHSS_change_from_d1_Day90[s1p.combined$patient == "s1p8"] <- NA  


# Raw absolute perihematomal edema volume (mL)
s1p.combined$PHE_Day1 <- NA
s1p.combined$PHE_Day8 <- NA
s1p.combined$PHE_Day15 <- NA

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p2"] <- 8.94
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p2"] <- 16.95
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p2"] <- 6.07

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p3"] <- 20.83
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p3"] <- 44.5
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p3"] <- 54.05

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p4"] <- 29.94
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p4"] <- 70.17
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p4"] <- 94.33

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p5"] <- 39.16
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p5"] <- 46.94
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p5"] <- 38.44

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p6"] <- 12.93
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p6"] <- 35.48
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p6"] <- 61.96

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p7"] <- 138.68 
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p7"] <- 170.75
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p7"] <- 176.85

s1p.combined$PHE_Day1[s1p.combined$patient == "s1p8"] <- 38.66
s1p.combined$PHE_Day8[s1p.combined$patient == "s1p8"] <- 33.12
s1p.combined$PHE_Day15[s1p.combined$patient == "s1p8"] <- 38.61

# Binary absolute perihematomal edema volume
s1p.combined$PHE_Day1_split_at_30 <- NA
s1p.combined$PHE_Day8_split_at_40 <- NA
s1p.combined$PHE_Day15_split_at_40 <- NA

s1p.combined$PHE_Day1_split_at_30[s1p.combined$PHE_Day1 >= 30] <- "PHE_Day1_above_or_equal_to_30"
s1p.combined$PHE_Day1_split_at_30[s1p.combined$PHE_Day1 < 30] <- "PHE_Day1_below_30"
s1p.combined$PHE_Day8_split_at_40[s1p.combined$PHE_Day8 >= 40] <- "PHE_Day8_above_or_equal_to_40"
s1p.combined$PHE_Day8_split_at_40[s1p.combined$PHE_Day8 < 40] <- "PHE_Day8_below_40"
s1p.combined$PHE_Day15_split_at_40[s1p.combined$PHE_Day15 >= 40] <- "PHE_Day15_above_or_equal_to_40"
s1p.combined$PHE_Day15_split_at_40[s1p.combined$PHE_Day15 < 40] <- "PHE_Day15_below_40"


# Raw hematoma volume (mL)
s1p.combined$Hematoma_Vol_Day1 <- NA
s1p.combined$Hematoma_Vol_Day8 <- NA
s1p.combined$Hematoma_Vol_Day15 <- NA

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p2"] <- 7.28
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p2"] <- 9.48
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p2"] <- 1.9

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p3"] <- 21.36
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p3"] <- 21.9
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p3"] <- 13.63

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p4"] <- 45.9
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p4"] <- 42.04
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p4"] <- 26.5

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p5"] <- 34.32
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p5"] <- 4.28
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p5"] <- 2.52

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p6"] <- 24.48
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p6"] <- 20.07
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p6"] <- 9.38

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p7"] <- 48.47
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p7"] <- 54.01
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p7"] <- 30.28

s1p.combined$Hematoma_Vol_Day1[s1p.combined$patient == "s1p8"] <- 26.63
s1p.combined$Hematoma_Vol_Day8[s1p.combined$patient == "s1p8"] <- 16.67
s1p.combined$Hematoma_Vol_Day15[s1p.combined$patient == "s1p8"] <- 5.56

# Binary hematoma volume outcome
s1p.combined$Hematoma_Vol_Day1_split_at_30 <- NA
s1p.combined$Hematoma_Vol_Day8_split_at_30 <- NA
s1p.combined$Hematoma_Vol_Day15_split_at_20 <- NA

s1p.combined$Hematoma_Vol_Day1_split_at_30[s1p.combined$Hematoma_Vol_Day1 >= 30] <- "Hematoma_Vol_Day1_above_or_equal_to_30"
s1p.combined$Hematoma_Vol_Day1_split_at_30[s1p.combined$Hematoma_Vol_Day1 < 30] <- "Hematoma_Vol_Day1_below_30"
s1p.combined$Hematoma_Vol_Day8_split_at_30[s1p.combined$Hematoma_Vol_Day8 >= 30] <- "Hematoma_Vol_Day8_above_or_equal_to_30"
s1p.combined$Hematoma_Vol_Day8_split_at_30[s1p.combined$Hematoma_Vol_Day8 < 30] <- "Hematoma_Vol_Day8_below_30"
s1p.combined$Hematoma_Vol_Day15_split_at_20[s1p.combined$Hematoma_Vol_Day15 >= 20] <- "Hematoma_Vol_Day15_above_or_equal_to_20"
s1p.combined$Hematoma_Vol_Day15_split_at_20[s1p.combined$Hematoma_Vol_Day15 < 20] <- "Hematoma_Vol_Day15_below_20"


# Raw relative perihematomal edema volume (mL)
s1p.combined$Rel_PHE_Day1 <- NA
s1p.combined$Rel_PHE_Day8 <- NA
s1p.combined$Rel_PHE_Day15 <- NA


# rPHE calculated based on ICH volume
s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p2"] <- 1.23
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p2"] <- 1.79
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p2"] <- 3.19

s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p3"] <- 0.974
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p3"] <- 2.03
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p3"] <- 3.97

s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p4"] <- 0.652
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p4"] <- 1.67
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p4"] <- 3.56

s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p5"] <- 1.14
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p5"] <- 10.97
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p5"] <- 15.25

s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p6"] <- 0.528
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p6"] <- 1.77
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p6"] <- 6.61

s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p7"] <- 2.86
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p7"] <- 3.16
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p7"] <- 5.84

s1p.combined$Rel_PHE_Day1[s1p.combined$patient == "s1p8"] <- 1.45
s1p.combined$Rel_PHE_Day8[s1p.combined$patient == "s1p8"] <- 1.99
s1p.combined$Rel_PHE_Day15[s1p.combined$patient == "s1p8"] <- 6.94

# Binary rPHE based on ICH volume
s1p.combined$Rel_PHE_Day1_split_at_1 <- NA
s1p.combined$Rel_PHE_Day8_split_at_3 <- NA
s1p.combined$Rel_PHE_Day15_split_at_5 <- NA

s1p.combined$Rel_PHE_Day1_split_at_1[s1p.combined$Rel_PHE_Day1 >= 1] <- "Rel_PHE_Day1_above_or_equal_to_1"
s1p.combined$Rel_PHE_Day1_split_at_1[s1p.combined$Rel_PHE_Day1 < 1] <- "Rel_PHE_Day1_below_1"
s1p.combined$Rel_PHE_Day8_split_at_3[s1p.combined$Rel_PHE_Day8 >= 3] <- "Rel_PHE_Day8_above_or_equal_to_3"
s1p.combined$Rel_PHE_Day8_split_at_3[s1p.combined$Rel_PHE_Day8 < 3] <- "Rel_PHE_Day8_below_3"
s1p.combined$Rel_PHE_Day15_split_at_5[s1p.combined$Rel_PHE_Day15 >= 5] <- "Rel_PHE_Day15_above_or_equal_to_5"
s1p.combined$Rel_PHE_Day15_split_at_5[s1p.combined$Rel_PHE_Day15 < 5] <- "Rel_PHE_Day15_below_5"


# Hemorrhage location
s1p.combined$Hematoma_Location <- NA

s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p2"] <- "Basal_ganglia"
s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p3"] <- "Basal_ganglia"
s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p4"] <- "Basal_ganglia"
s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p5"] <- "Occipital_lobe"
s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p6"] <- "Basal_ganglia"
s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p7"] <- "Occipital_lobe" 
s1p.combined$Hematoma_Location[s1p.combined$patient == "s1p8"] <- "Occipital_lobe"


# ICH Score
s1p.combined$ICH_score <- NA

s1p.combined$ICH_score[s1p.combined$patient == 's1p2'] <- 0
s1p.combined$ICH_score[s1p.combined$patient == 's1p3'] <- 0
s1p.combined$ICH_score[s1p.combined$patient == 's1p4'] <- 2
s1p.combined$ICH_score[s1p.combined$patient == 's1p5'] <- 1
s1p.combined$ICH_score[s1p.combined$patient == 's1p6'] <- 1
s1p.combined$ICH_score[s1p.combined$patient == 's1p7'] <- 2
s1p.combined$ICH_score[s1p.combined$patient == 's1p8'] <- 0
  
# Binary ICH score
s1p.combined$ICH_score_split_at_1 <- NA
s1p.combined$ICH_score_split_at_1[s1p.combined$ICH_score >= 1] <- "ICH_score_above_or_equal_to_1"
s1p.combined$ICH_score_split_at_1[s1p.combined$ICH_score < 1] <- "ICH_score_below_1"

s1p.combined$ICH_score_split_at_2 <- NA
s1p.combined$ICH_score_split_at_2[s1p.combined$ICH_score >= 2] <- "ICH_score_above_or_equal_to_2"
s1p.combined$ICH_score_split_at_2[s1p.combined$ICH_score < 2] <- "ICH_score_below_2"


save(s1p.combined, file = "s1p.combined_metadata.RData")

