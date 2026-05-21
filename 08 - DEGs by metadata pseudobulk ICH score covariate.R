library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
if (!require(Matrix.utils)) remotes::install_github("cvarrichio/Matrix.utils"); library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)
library(limma)
library(purrr)
library(EnhancedVolcano)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)


mywd <- paste0("[path]", myTime(), " DEGs by metadata correct_HV ICH_score covariate removed corrected contrasts mincount=10 mintotalcount=100")
dir.create(mywd)
setwd(mywd)


#Load Clustered Seurat Object
load("[path]/s1p.combined_correct_HV.RData")

s1p.combined <- subset(s1p.combined, lineageClusterIdents == 'CD4+ T cells' |
                         lineageClusterIdents == 'CD8+ T cells' |
                         lineageClusterIdents == 'B cells' |
                         lineageClusterIdents == 'NK cells' |
                         lineageClusterIdents == 'Monocytes' |
                         lineageClusterIdents == 'CD16+ Monocytes')

# Pseudobulk pipeline based on https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pseudobulkDE_edgeR.md

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = s1p.combined, layer = "counts", assay = "RNA")
metadata <- s1p.combined@meta.data
Idents(s1p.combined) <- "lineageClusterIdents"
metadata$cluster_id <- metadata$lineageClusterIdents
metadata$sample_id <- paste(metadata$patient, metadata$time)
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)


## Remove lowly expressed genes which have less than 10 cells with any counts
dim(sce)
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Initialize lists to store DEGs
metadataEdgeRDEGs <- list()
metadataEdgeRDEGTStatRanks <- list()
metadataEdgeRDEGTStatRanksVector <- list()


# Function to calculate DEGs for a given metadata field at a given timepoint
calcDEGs <- function(sceDEG, metadataField, metadataName, day, myContrast) {
  
  # Named vector of cluster names
  kids <- purrr::set_names(levels(as.factor(sceDEG$lineageClusterIdents)))
  nk <- length(kids)
  sids <- purrr::set_names(levels(as.factor(sceDEG$sample_id)))
  ns <- length(sids)

  print(metadataName)
  print(day)
  print(myContrast)
  class(myContrast)
  # Generate sample level metadata
  table(sceDEG$sample_id)
  n_cells <- table(sceDEG$sample_id) %>%  as.vector()
  names(n_cells) <- names(table(sceDEG$sample_id))
  m <- match(names(n_cells), sceDEG$sample_id)
  sceDEG$group <- metadataField # Define the comparison metadata
  ei <- data.frame(colData(sceDEG)[m, ], 
                   n_cells, row.names = NULL) 
  kable(ei)
  
  # Aggregate the counts per sample_id and cluster_id to make the pseudobulk object
  groups <- colData(sceDEG)[, c("cluster_id", "sample_id")]
  groups$sample_id <- factor(groups$sample_id)
  pb <- aggregate.Matrix(t(counts(sceDEG)), 
                         groupings = groups, fun = "sum") 
  pb[1:8, 1:8]
  splitf <- sapply(stringr::str_split(rownames(pb), 
                                      pattern = "_",n = 2), `[`, 1)
  pb <- split.data.frame(pb,factor(splitf)) %>%
    lapply(function(u) 
      set_colnames(t(u), gsub(".*_", "", rownames(u))))
  
  # Print out the table of cells in each cluster-sample group
  options(width = 100)
  kable(table(sceDEG$cluster_id, sceDEG$sample_id))
  
  # Select the clusters to include
  keepClusters <- as.character(unique(sceDEG$cluster_id))
  keepClusters <- keepClusters[keepClusters %in% c('CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells', 'Monocytes', 'CD16+ Monocytes')]
  (interestingClusters <- SingleCellExperiment(assays = pb[keepClusters]))

  # compute MDS (similar to PCA)
  mds <-  lapply(as.list(assays(interestingClusters)), function(a){
    DGEList(a, remove.zeros = TRUE) %>% 
      calcNormFactors %>% 
      plotMDS.DGEList(plot = FALSE)
  })
  cnames <- paste("Cluster", keepClusters) 
  for (m in 1:length(mds)){
    mds[[m]]$cluster <- cnames[m]
  }
  plots <- lapply(mds, function(m){
    gg_df <- data.frame(m[c("x", "y")],
                        sample_id = sids,
                        group_id = ei$group,
                        cluster_id = rep(m$cluster, length(m$x)))})
  plotFunc <- function(x) {
    ggplot(x, aes(x, y, col = group_id)) + 
      geom_point(size = 3, alpha = 0.8) +
      labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
      ggtitle(unique(x$cluster_id)) + 
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      coord_fixed() 
  }
  do.call(grid.arrange,c(lapply(plots, plotFunc)))
  
  # Replot mds but faceting within ggplot
  plots_df <- as.data.frame(plots[[1]])
  colnames(plots_df) <- colnames(plots[[1]])
  for (i in 2:length(plots)) {
    plots_df <- rbind(plots_df, plots[[i]])
  }
  
  plots_df$cluster_id <- as.factor(plots_df$cluster_id)
  
  ggplot(plots_df, aes(x, y, color = group_id, group = cluster_id)) + 
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
    ggtitle(paste0('MDS by ', metadataName, ' ', day)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ cluster_id)
  ggsave(paste0('MDS by ', metadataName, ' ', day, ' ICH_score covariate removed ', myTime(), '.pdf'), height = 5, width = 7)
  
  
  #Voom LM Fit Model 
  (design <- model.matrix(~ 0 + ei$group + ei$ICH_score) %>% 
      set_rownames(ei$sample_id) %>% 
      set_colnames(c(levels(factor(ei$group)), 'ICH_score')))
  
  (contrast <- makeContrasts(contrasts = myContrast, levels = colnames(design)))
  
  
  # Use voom to create list of genes for GSEA, ranked by t-statistic
  dflist <- list()
  fit <- list()
  efit <- list()
  lcpm <- list()
  # res <- list()
  for (i in 1:length(keepClusters)) {
    dflist[[i]] = DGEList(assays(interestingClusters)[[i]])
    dflist[[i]] = calcNormFactors(dflist[[i]], method = "RLE")
    genes.use1 = filterByExpr(dflist[[i]], design = design, min.count = 10, min.total.count = 100)
    dflist[[i]] = dflist[[i]][ genes.use1, keep.lib.sizes=FALSE]
    fit[[i]] <- voomLmFit(dflist[[i]],  design, plot=FALSE, sample.weights=TRUE)
    fit[[i]] <- contrasts.fit(fit[[i]], 
                              #coefficients = colnames(design)[1],
                              contrasts = contrast)
    efit[[i]] <- eBayes(fit[[i]])
    plotSA(efit[[i]])
  }
  
  ranks = list()
  for (x in 1:length(efit)) {
    test = as.data.frame(topTable(efit[[x]], coef = NULL, number = Inf))
    ranks[[x]] =
      test %>% rownames_to_column("gene") %>%
      arrange(desc(t)) %>%
      select(gene, t) %>%
      column_to_rownames("gene") %>%
      t() %>%
      unlist(use.names = T)
    ranks[[x]] = ranks[[x]][1, ]
  }
  
  # DGE Analysis with EdgeR
  res <- lapply(keepClusters, function(k) {
    y <- assays(interestingClusters)[[k]]
    y <- DGEList(y, remove.zeros = TRUE)
    keep <- filterByExpr(y, design = design, group = NULL, lib.size = NULL, min.count = 10, min.total.count = 100)
    y <- y[keep,]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    fit <- glmQLFTest(fit, contrast = contrast)
    topTags(fit, n = Inf, sort.by = "none")$table %>%
      dplyr::mutate(gene = rownames(.), cluster_id = k) %>%
      dplyr::rename(p_val = PValue, p_adj = FDR)
  })
  
  # Results filtering & overview: filter FDR < 0.05, |logFC| > 1 & sort by FDR
  res_fil <- lapply(res, 
                    function(u)  u %>% 
                      dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                      dplyr::arrange(p_adj))
  
  ## Count the number of differential findings by cluster.
  n_de <- vapply(res_fil, nrow, numeric(1))
  cbind(cluster=keepClusters, numDE_genes=n_de, 
        percentage = round(n_de / nrow(interestingClusters) * 100, digits =2)) %>%  kable()
  
  # Output results
  for(i in 1:length(keepClusters)){
    cluster <- unique(res[[i]]$cluster_id)
    
    #Write DEGs as CSVs
    out <- res[[i]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
    write.csv(out, file = paste0("EdgeR DEGs by ", metadataName, " ", myContrast, " ", cluster, " ", day, " ICH_score covariate removed ", myTime(), ".csv"), quote=F, row.names = F)
    metadataEdgeRDEGs[[paste(metadataName, myContrast, cluster, day)]] <- out
    
    # Write ranked list as .txt and .rnk (for GSEA)
    out <- as.table(ranks[[i]])
    write.table(out, file = paste0("EdgeR TStatRanks by ", metadataName, " ", myContrast, " ", cluster, " ", day, " ICH_score covariate removed ", myTime(), ".txt"), 
                quote=F, row.names = F, col.names = F, sep = "\t")

    # Save ranked list as data.frame
    out <- as.data.frame(out)
    colnames(out) <- c("Gene", "TStatRank")
    metadataEdgeRDEGTStatRanks[[paste(metadataName, myContrast, cluster, day)]] <- out
    
    # Save ranked list as named numeric vector
    outVector <- setNames(out$TStatRank, out$Gene)
    metadataEdgeRDEGTStatRanksVector[[paste(metadataName, myContrast, cluster, day)]] <- outVector
  }
  gc()
  return(list(metadataEdgeRDEGs, metadataEdgeRDEGTStatRanks, metadataEdgeRDEGTStatRanksVector))
}


for (d in c("Day 1", "Day 3", "Day 7")) {
  sceTemp <- sce[,sce$time == d]

  tempList <- calcDEGs(sceTemp, sceTemp$MRS_Day90_split_at_4, "MRS_Day90_split_at_4", d, "MRS_Day90_below_4 - MRS_Day90_above_or_equal_to_4")
  metadataEdgeRDEGs <- tempList[[1]]
  metadataEdgeRDEGTStatRanks <- tempList[[2]]
  metadataEdgeRDEGTStatRanksVector <- tempList[[3]]
}

save(metadataEdgeRDEGs, metadataEdgeRDEGTStatRanks, metadataEdgeRDEGTStatRanksVector, file = "EdgeR_Metadata_DEG_Lists_covariate_ICH_score.RData")

save.image(file = "EdgeR_Metadata_DEG_Lists_covariate_ICH_score_Seurat.RData")



