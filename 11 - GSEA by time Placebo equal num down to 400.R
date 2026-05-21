library(dplyr)
library(fgsea)
library(ggplot2)
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

mywd <- paste0("[path]/", myTime(), " scDEG by time Placebo equal num down to 400 cells")
dir.create(mywd)
setwd(mywd)

# (Load datasets below, directly before running GSEA)

# Code from https://biostatsquid.com/fgsea-tutorial-gsea/
# Function: Adjacency matrix to list
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Graph all three timepoint comparisons
graphBarChartHorizontal <- function(myData, cellType) {
  if (nrow(myData) > 0) { # Make sure there's data in the dataframe
    
    myData <- filter(myData, str_detect(pathway, paste(paste0('^', myPathways, '$'), collapse = "|")))
    
    myData$pathway <- factor(myData$pathway, levels = rev(myPathways))
    
    myData$timepoint <- factor(myData$timepoint, levels = c("D3vD1", "D7vD1", "D7vD3"))
    
    # Find width of longest canonical pathway label
    cpLabels <- as.character(myData$pathway)
    longestCPLabel <- max(nchar(cpLabels))
    
    colors = setNames(hue_pal()(3), c("D3vD1", "D7vD1", "D7vD3"))
    
    ggplot(myData, aes(x = pathway, y = NES, fill = timepoint, color = timepoint)) + 
      geom_bar(aes(alpha = -log10(padj)), stat = "identity", linewidth = 0.4, width = 0.8, position = position_dodge(width = 0.9)) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = .5),
            axis.title = element_text(size = 15)) +
      labs(y = 'NES', x = 'GSEA Hallmark Pathways',
           alpha = '-log(adj p-value)', color = 'Timepoint', fill = 'Timepoint',
           title = paste(cellType, "- Top pathways (scDEGs)\nscDEGs equalNum\nDown to 400 cells\nmin.pct=0.1\n")) +
      scale_y_continuous(limits = c(ifelse(min(myData$NES) < 0, min(myData$NES)*1.5, 0), ifelse(max(myData$NES) > 0, max(myData$NES)*1.5, 0)))
    ggsave(paste0("Barchart Horizontal GSEA Hallmark scDEGs by time Placebo equalNum down to 400 cells min_pct=0_1 ", cellType, " " , myTime(),  ".pdf"), 
           width = (3.3 + .35 * length(myPathways)), height = (4 + 0.06 * longestCPLabel), device = cairo_pdf) 
  }
}

# Graph D3vD1 and D7vD1
graphBarChartHorizontal_D3vD1_D7vD1 <- function(myData, cellType) {
  if (nrow(myData) > 0) { # Make sure there's data in the dataframe
    
    myData <- filter(myData, str_detect(pathway, paste(paste0('^', myPathways, '$'), collapse = "|")))
    
    # Plot only D3vD1 and D7vD1
    myData <- filter(myData, timepoint != "D7vD3") 
   
    myData$pathway <- factor(myData$pathway, levels = rev(myPathways))
    
    myData$timepoint <- factor(myData$timepoint, levels = c("D3vD1", "D7vD1", "D7vD3"))
    
    # Find width of longest canonical pathway label
    cpLabels <- as.character(myData$pathway)
    longestCPLabel <- max(nchar(cpLabels))
    
    colors = setNames(hue_pal()(3), c("D3vD1", "D7vD1", "D7vD3"))
    
    ggplot(myData, aes(x = pathway, y = NES, fill = timepoint, color = timepoint)) + 
      geom_bar(aes(alpha = -log10(padj)), stat = "identity", linewidth = 0.4, width = 0.8, position = position_dodge(width = 0.9)) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      theme_classic() + 
      theme(axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = .5),
            axis.title = element_text(size = 15)) +
      labs(y = 'NES', x = 'GSEA Hallmark Pathways',
           alpha = '-log(adj p-value)', color = 'Timepoint', fill = 'Timepoint',
           title = paste(cellType, "- Top pathways (scDEGs)\nscDEGs equalNum\nDown to 400 cells\nmin.pct=0.1\n")) +
      scale_y_continuous(limits = c(ifelse(min(myData$NES) < 0, min(myData$NES)*1.5, 0), ifelse(max(myData$NES) > 0, max(myData$NES)*1.5, 0)))
    ggsave(paste0("Barchart Horizontal GSEA Hallmark scDEGs by time Placebo D3vD1 and D7vD1 equalNum down to 400 cells min_pct=0_1 ", cellType, " " , myTime(),  ".pdf"), 
           width = (3.3 + .25 * length(myPathways)), height = (4 + 0.06 * longestCPLabel), device = cairo_pdf) 
  }
}


# Function to perform GSEA (Code from https://bioinformaticsbreakdown.com/how-to-gsea/ (Brian Gudenas))
GSEA = function(gene_list, GMT_file, pval, myGeneset, myComparison, myTimepoint, myCluster) {
  set.seed(1111)

  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGMT = fgsea::gmtPathways(GMT_file)
  
  # Subset GMT lists on only the genes present in our dataset to avoid bias (from https://biostatsquid.com/fgsea-tutorial-gsea/)
  GMTgenes <- unique(unlist(myGMT))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  GMTmat <- matrix(NA, dimnames = list(GMTgenes, names(myGMT)),
                nrow = length(GMTgenes), ncol = length(myGMT))
  for (i in 1:dim(GMTmat)[2]){
    GMTmat[,i] <- as.numeric(GMTgenes %in% myGMT[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  intersectGenes <- intersect(names(gene_list), GMTgenes)
  GMTmat <- GMTmat[intersectGenes, colnames(GMTmat)[which(colSums(GMTmat[intersectGenes,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again using the function we previously defined
  final_list <- matrix_to_list(GMTmat)
  myGMT <- final_list
  
  fgRes <- fgsea::fgsea(pathways = myGMT,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=400#, ## maximum gene set size
                        ) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGMT,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  
  # Clean up pathway names (specific for Hallmark pathways)
  fgRes$pathway <- str_sub(fgRes$pathway, 10) %>%
    str_replace_all("_", " ")
  
  # Filter on pathways with highest absolute NES
  orderedRows <- order(abs(fgRes$NES), decreasing = TRUE)
  orderedRows_top10 <- orderedRows[1:10][!is.na(orderedRows[1:10])]
  filtRes10 = fgRes[orderedRows_top10,]
  
  # Placeholder plot in case there are no pathways to graph
  g10 <- ggplot()
  
  if (nrow(filtRes10) > 0) {
    total_up = sum(fgRes$Enrichment == "Up-regulated")
    total_down = sum(fgRes$Enrichment == "Down-regulated")
    header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
    
    colors = setNames(c("red", "blue"),
                     c("Up-regulated", "Down-regulated"))
    
    # find x min and max to manually set axis limits
    xmin <- min(filtRes10$NES)
    xmax <- max(filtRes10$NES)
    myRange <- xmax - xmin
    
    # Find width of longest canonical pathway label
    myLabels <- as.character(filtRes10$pathway)
    longestLabel <- max(nchar(myLabels))
    
    g10= 
      ggplot(filtRes10, aes(NES, reorder(pathway, NES))) +
        theme_classic() +
        geom_point(aes(color = Enrichment, size = -log10(padj)), alpha = 0.6) +
        scale_color_manual(values = colors) +
        scale_size_continuous(range = c(2,8)) +
        theme(axis.text.y = element_text(size = 12),
              axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
              axis.title = element_text(size = 15),
              title = element_text(size = 10)) +
        labs(y = paste('GSEA', myGeneset), x = 'NES',
             size = '-log(FDR q-val)', color = 'Direction',
             title = paste0("GSEA NoNPerm top 10\nPlacebo by Time - scDEGs\nscDEGs equalNum; down to 400 cells; min.pct=0.1\n", myCluster, " ", myTimepoint, "\n", myGeneset, "\n ")) +
        xlim(xmin - 0.2 * myRange, xmax + 0.2 * myRange)
      
    g10
    ggsave(paste0("Dotplot GSEA scDEGs NoNPerm top 10 Placebo by Time equalNum down to 400 cells min_pct=0_1 ", myGeneset, " ", myCluster, " ", myTimepoint, " ", myTime(),  ".pdf"), 
           height = (2.2 + .25 * nrow(filtRes10)), width = (3.2 + 0.11 * longestLabel), device = cairo_pdf)
  }
  
  # Same, but plot top 20
  orderedRows_top20 <- orderedRows[1:20][!is.na(orderedRows[1:20])]
  filtRes20 = fgRes[orderedRows_top20,]
  
  # Placeholder plot in case there are no pathways to graph
  g20 <- ggplot()
  
  if (nrow(filtRes20) > 0) {
    total_up = sum(fgRes$Enrichment == "Up-regulated")
    total_down = sum(fgRes$Enrichment == "Down-regulated")
    header = paste0("Top 20 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
    
    colors = setNames(c("red", "blue"),
                      c("Up-regulated", "Down-regulated"))
    
    # find x min and max to manually set axis limits
    xmin <- min(filtRes20$NES)
    xmax <- max(filtRes20$NES)
    myRange <- xmax - xmin
    
    # Find width of longest canonical pathway label
    myLabels <- as.character(filtRes20$pathway)
    longestLabel <- max(nchar(myLabels))
    
    # Make dot plot
    g20= 
      ggplot(filtRes20, aes(NES, reorder(pathway, NES))) +
      theme_classic() +
      geom_point(aes(color = Enrichment, size = -log10(padj)), alpha = 0.6) +
      scale_color_manual(values = colors) +
      scale_size_continuous(range = c(2,8)) +
      theme(axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 0.95),
            axis.title = element_text(size = 15),
            title = element_text(size = 10)) +
      #geom_vline(xintercept = 0) +
      labs(y = paste('GSEA', myGeneset), x = 'NES',
           size = '-log(FDR q-val)', color = 'Direction',
           title = paste0("GSEA NoNPerm top 20\nPlacebo by Time - scDEGs\nscDEGs equalNum; down to 400 cells; min.pct=0.1\n", myCluster, " ", myTimepoint, "\n", myGeneset, "\n ")) +
      xlim(xmin - 0.2 * myRange, xmax + 0.2 * myRange)
    
    g20
    ggsave(paste0("Dotplot GSEA scDEGs NoNPerm top 20 Placebo by Time equalNum down to 400 cells min_pct=0_1 ", myGeneset, " ", myCluster, " ", myTimepoint, " ", myTime(),  ".pdf"), 
           height = (2.2 + .25 * nrow(filtRes20)), width = (3.2 + 0.11 * longestLabel), device = cairo_pdf)
  }
  
  output = list("Results" = fgRes, "Plot_top10" = g10, "Plot_top20" = g20)
  return(output)
}

# Build list of scDEGs calculated with Seurat
scDEGListbyTimePlacebo <- list()
clusters <- c('Monocytes', 'CD16+ Monocytes', 'CD4+ T cells', 'CD8+ T cells', 'NK cells', 'B cells')
timepoints <- c('Day 3 vs Day 1', 'Day 7 vs Day 1', 'Day 7 vs Day 3')
for (c in clusters) {
  for (t in timepoints) {
    # Load DEGs calculated in script "08 - DEGs by time Placebo-treated equal num down to 400.R"
    filename <- paste0('[path]/DEGs Placebo-treated equalNum down to 400 cells ', c, ' ', t, '.txt')
    if (file.exists(filename)){
      df <- read_delim(filename, delim = '\t')
      df <- df[order(df$avg_log2FC, decreasing = TRUE),]
      vector <- df$avg_log2FC
      names(vector) <- df$...1
      t <- case_when(
        t == 'Day 3 vs Day 1' ~ 'D3vD1',
        t == 'Day 7 vs Day 1' ~ 'D7vD1',
        t == 'Day 7 vs Day 3' ~ 'D7vD3'
      )
      scDEGListbyTimePlacebo[[paste('time', c, t)]] <- vector
    }
  }
}

# Calculate and graph bubble plots
GSEA_res_Time_Hallmark <- list()
for (i in 1:length(scDEGListbyTimePlacebo)) {
  gene_list <- scDEGListbyTimePlacebo[[i]]
  GMT_file = "[path]/msigdb_v2023.2.Hs_GMTs/symbols/h.all.v2023.2.Hs.symbols.gmt"
  
  geneset <- "MSigDB Hallmark Geneset"
  dataset <- str_split(names(scDEGListbyTimePlacebo)[i], pattern = " ")[[1]]
  comparison <- dataset[1]
  timepoint <- str_c(dataset[length(dataset)], collapse = " ")
  dataset <- dataset[c(-1, -length(dataset))]
  cluster <- str_c(dataset, collapse = " ")
  
  res = GSEA(gene_list, GMT_file, pval = 0.05, geneset, comparison, timepoint, cluster)
  GSEA_res_Time_Hallmark[[names(scDEGListbyTimePlacebo)[i]]] <- res
  gc()
}
save(GSEA_res_Time_Hallmark, file = "GSEA_scDEGs_equalNum_down_to_400_cells_Time_Placebo_Hallmark.RData")

# Load previously run GSEA results
#load("[path]/GSEA_scDEGs_equalNum_down_to_400_cells_Time_Placebo_Hallmark.RData")

# Plot all timepoints
for (pop in c('CD4+ T cells', 'CD8+ T cells', 'NK cells', 'B cells', 'Monocytes', 'CD16+ Monocytes')) {
  d3vd1 <- GSEA_res_Time_Hallmark[[paste('time', pop, 'D3vD1')]]$Results
  d7vd1 <- GSEA_res_Time_Hallmark[[paste('time', pop, 'D7vD1')]]$Results
  d7vd3 <- GSEA_res_Time_Hallmark[[paste('time', pop, 'D7vD3')]]$Results
 
  # To allow empty dfs to be combined
  if (nrow(d3vd1) > 0) {
    d3vd1$timepoint <- "D3vD1"
  }
  if (nrow(d7vd1) > 0) {
    d7vd1$timepoint <- "D7vD1"
  }
  if (nrow(d7vd3) > 0) {
    d7vd3$timepoint <- "D7vD3"
  }
  
  class(d3vd1$Enrichment) <- 'character'
  class(d7vd1$Enrichment) <- 'character'
  class(d7vd3$Enrichment) <- 'character'
  
  barchartData <- bind_rows(d3vd1, d7vd1) %>%
    bind_rows(d7vd3)
  
  # Fill in missing pathways
  bcDataSpread <- select(barchartData, c(pathway, timepoint, NES))
  bcDataSpread <- spread(bcDataSpread, key = timepoint, value = NES)
  if (!("D3vD1" %in% colnames(bcDataSpread))) { # Add timepoints that have no pathways at all
    bcDataSpread$D3vD1 <- 0
  }
  if (!("D7vD1" %in% colnames(bcDataSpread))) {
    bcDataSpread$D7vD1 <- 0
  }
  if (!("D7vD3" %in% colnames(bcDataSpread))) {
    bcDataSpread$D7vD3 <- 0
  }
  bcDataGather <- gather(bcDataSpread, key = timepoint, value = NES, 2:ncol(bcDataSpread))
  barchartData <- left_join(bcDataGather, barchartData, by = c('pathway', 'timepoint'), suffix = c("", ".y")) %>%
    select(-ends_with(".y")) # join dfs and remove duplicate column
  barchartData[is.na(barchartData)] <- 0 # replace NAs with zeros
  
  # Select top 15 pathways
  barchartData <- barchartData[order(abs(barchartData$NES), decreasing = TRUE),]
  barchartData <- filter(barchartData, pathway %in% unique(barchartData$pathway)[1:15])
  
  # Order by NES
  barchartData <- barchartData[order(barchartData$NES),]
  if (nrow(filter(barchartData, timepoint == 'D3vD1')) > 0) {
    D3vD1only <- filter(barchartData, timepoint == 'D3vD1')
    D3vD1order <- D3vD1only[order(D3vD1only$NES),]$pathway
    D3vD1order <- factor(D3vD1order, levels = D3vD1order)
    barchartData$pathway <- factor(barchartData$pathway, levels = D3vD1order)
  } else if (nrow(filter(barchartData, timepoint == 'D7vD1')) > 0) {
    D7vD1only <- filter(barchartData, timepoint == 'D7vD1')
    D7vD1order <- D7vD1only[order(D7vD1only$NES),]$pathway
    D7vD1order <- factor(D7vD1order, levels = D7vD1order)
    barchartData$pathway <- factor(barchartData$pathway, levels = D7vD1order)
  } else if (nrow(filter(barchartData, timepoint == 'D7vD3')) > 0) {
    D7vD3only <- filter(barchartData, timepoint == 'D7vD3')
    D7vD3order <- D7vD3only[order(D7vD3only$NES),]$pathway
    D7vD3order <- factor(D7vD3order, levels = D7vD3order)
    barchartData$pathway <- factor(barchartData$pathway, levels = D7vD3order)
  }
  
  myPathways <- unique(barchartData$pathway)
  graphBarChartHorizontal(barchartData, pop)
  
  # Plot D3vD1 and D7vD1 only
  if (nrow(filter(barchartData, timepoint == 'D3vD1')) > 0) {
    D3vD1only <- filter(barchartData, timepoint == 'D3vD1')
    D3vD1order <- D3vD1only[order(D3vD1only$NES),]$pathway
    D3vD1order <- factor(D3vD1order, levels = D3vD1order)
    barchartData$pathway <- factor(barchartData$pathway, levels = D3vD1order)
  } else if (nrow(filter(barchartData, timepoint == 'D7vD1')) > 0) {
    D7vD1only <- filter(barchartData, timepoint == 'D7vD1')
    D7vD1order <- D7vD1only[order(D7vD1only$NES),]$pathway
    D7vD1order <- factor(D7vD1order, levels = D7vD1order)
    barchartData$pathway <- factor(barchartData$pathway, levels = D7vD1order)
  }
  
  # Filter out pathways that were only significant in D7vD3 comparison
  barchartDataNoZeros <- filter(barchartData, timepoint != 'D7vD3')
  barchartDataNoZeros <- select(barchartDataNoZeros, c('pathway', 'timepoint', 'NES'))
  barchartDataNoZeros <- spread(barchartDataNoZeros, key = timepoint, value = NES)
  barchartDataNoZeros <- mutate(barchartDataNoZeros, 'mySum' = abs(D3vD1) + abs(D7vD1))
  barchartDataNoZeros <- filter(barchartDataNoZeros, mySum != 0)
  
  myPathways <- factor(barchartDataNoZeros$pathway)
  
  graphBarChartHorizontal_D3vD1_D7vD1(barchartData, pop)
}








