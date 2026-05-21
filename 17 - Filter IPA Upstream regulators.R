library(tidyverse)
library(stringr)

# Set wd to the path containing IPA upstream regulator output from Ingenuity Pathway Analysis
setwd("[path]")

file_list <- list.files(path = getwd())
file_list <- file_list[grep("txt$", file_list)]

# Remove Qiagen copyright
for (i in 1:length(file_list)) {
  temp_data <- read.delim(file_list[i], header = FALSE, sep = "\t", check.names = FALSE)
  temp_data <- temp_data[-1,]
  colnames(temp_data) <- temp_data[1,]
  temp_data <- temp_data[-1,]
  write.table(temp_data, file = paste0("./No copyright/", file_list[i]), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

setwd("[path]/No copyright")

file_list <- list.files(path = getwd())

# Filter on p <0.05 and number of molecules >4
for (i in 1:length(file_list)) {
  temp_data <- read.delim(file_list[i], header = TRUE, sep = "\t", check.names = FALSE)
  header <- colnames(temp_data)
  temp_data <- read.delim(file_list[i], header = TRUE, sep = "\t")
  temp_data <- filter(temp_data, B.H.corrected.p.value < 0.05)
  colnames(temp_data) <- header
  
  # Filter on number of molecules
  temp_data$Molecule_count <- ""
  fun_count <- function(x) {
    return (str_count(x, pattern = ",") + 1)
  }
  temp_data <- mutate(temp_data, Molecule_count = fun_count(`Target Molecules in Dataset`))
  temp_data <- filter(temp_data, Molecule_count > 4)
  
  write.table(temp_data, file = paste0("../Filtered p<0_05 >4 molecules/", substr(file_list[i], 1, nchar(file_list[i]) - 4), " ", nrow(temp_data), " upstream regulators.txt"), 
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

# Filter on cytokines
for (i in 1:length(file_list)) {
  # Filter on p-value
  temp_data <- read.delim(file_list[i], header = TRUE, sep = "\t", check.names = FALSE)
  header <- colnames(temp_data)
  temp_data <- read.delim(file_list[i], header = TRUE, sep = "\t")
  temp_data <- filter(temp_data, B.H.corrected.p.value < 0.05)
  colnames(temp_data) <- header
  
  # Filter on number of molecules
  if (nrow(temp_data) > 0) {
    temp_data$Molecule_count <- ""
    fun_count <- function(x) {
      return (str_count(x, pattern = ",") + 1)
    }
    
    temp_data <- mutate(temp_data, Molecule_count = fun_count(`Target Molecules in Dataset`))
    temp_data <- filter(temp_data, Molecule_count > 4)
  }
  
  # Filter on cytokines
  temp_data <- filter(temp_data, `Molecule Type` == "cytokine")
  colnames(temp_data) <- header
  
  write.table(temp_data, file = paste0("../Filtered p<0_05 >4 molecules Cytokines/", substr(file_list[i], 1, nchar(file_list[i]) - 4), " ", nrow(temp_data), " upstream regulators.txt"), 
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}




