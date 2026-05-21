library(tidyverse)

# Set wd to the path of the DEGs calculated in script "08 - DEGs by treatment equal num down to 1600.R"
setwd("[path]")

file_list <- list.files(path = getwd())
file_list <- file_list[grep("txt$", file_list)]

for (i in 1:length(file_list)) {
  temp_data <- read.table(file_list[i], header = TRUE, sep = "\t")
  temp_data <- select(temp_data, X, avg_log2FC, p_val_adj)
  colnames(temp_data) <- c("gene", "avg_log2FC", "p_value_adj")
  write.table(temp_data, file = paste("./For IPA/For IPA", file_list[i]), 
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}
