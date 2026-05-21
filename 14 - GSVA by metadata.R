library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(GSEABase)
library(GSVA)
library(viridis)
library(msigdb)
library(ExperimentHub)

# Function to return current date and time in a filename-friendly format
myTime <- function() {
  t <- gsub(":", "", Sys.time())
  t <- gsub(" ", "_", t)
  t <- gsub("\\.\\d+", "", t)
  return(t)
}

set.seed(1111)

mywd <- paste0("[path]/", myTime(), ' by metadata')
dir.create(mywd)
setwd(mywd)

# Load clustered data with metadata added from Script 03 - Add Patient Metadata
load("[path]/s1p.combined_metadata.RData")

DefaultAssay(s1p.combined) <- "RNA"

s1p.combined <- NormalizeData(s1p.combined)

s1p.combined <- SetIdent(s1p.combined, value = "lineageClusterIdents")

s1p.subset <- s1p.combined

# Rename BAF312_interrupted to remove the underscore
s1p.subset$treatment <- str_replace_all(s1p.subset$treatment, "_", "-")

# Save metadata to add to pb object later
meta <- s1p.subset@meta.data %>%
  unique() %>%
  remove_rownames()
meta$time.treatment.patient.lineage <- str_c(meta$time, meta$treatment, meta$patient, meta$lineageClusterIdents, sep = "_")
View(dplyr::select(meta, c('time', 'treatment', 'patient', 'lineageClusterIdents', 'time.treatment.patient.lineage')))
meta <- meta[,25:70]

# Create pseudobulk object
pb <- AggregateExpression(s1p.subset, assays = "RNA", return.seurat = TRUE, 
                          group.by = c("time", "treatment", "patient", "lineageClusterIdents"))
pb$time.treatment.patient.lineage <- names(Idents(pb))

# Create metadata fields for pb
pb$time <- NA
pb$treatment <- NA
pb$patient <- NA
pb$lineageClusterIdents <- NA
t <- str_split(names(pb$time.treatment.patient.lineage), "_")
for (i in 1:ncol(pb)) {
  pb$time[i] <- t[[i]][1]
  pb$treatment[i] <- t[[i]][2]
  pb$patient[i] <- t[[i]][3]
  pb$lineageClusterIdents[i] <- t[[i]][4]
}

# Number of cells by sample and celltype
n_cells <- s1p.subset@meta.data %>% 
  dplyr::count(time, treatment, patient, lineageClusterIdents) %>%
  rename("n" = "n_cells")

meta_pb <- left_join(pb@meta.data, n_cells)
rownames(meta_pb) <- meta_pb$time.treatment.patient.lineage
pb@meta.data <- meta_pb

meta <- meta[!duplicated(meta),]
pb@meta.data <- left_join(pb@meta.data, meta, by = "time.treatment.patient.lineage")
rownames(pb@meta.data) <- pb$time.treatment.patient.lineage


colnames(GetAssayData(pb, assay = "RNA", slot = "counts"))
rownames(pb@meta.data)

pb.subset <- subset(pb, subset = lineageClusterIdents == "Monocytes")

#save.image('Pseudobulk_for_GSEA.RData')
load('[path]/Pseudobulk_for_GSEA.RData')



# Build MSigDB geneset lists to plot together, per https://bioconductor.org/packages//release/data/experiment/vignettes/msigdb/inst/doc/msigdb.html#download-data-from-the-msigdb-r-package
eh = ExperimentHub()
query(eh, 'msigdb')

msigdb.hs <- getMsigdb(org = 'hs', id = 'SYM')
names(msigdb.hs)

# Hallmark Genesets
gs.hallmark <- msigdb.hs[grep("^HALLMARK", names(msigdb.hs))]
gs.hallmark.list <- list()
for (i in 1:length(gs.hallmark)) {
  gs.hallmark.list[[names(gs.hallmark[i])]] <- gs.hallmark[[i]]@geneIds
}
gs <- gs.hallmark

# Run GSVA
normalizedcounts <- pb.subset@assays$RNA$data

gsvaPar <- gsvaParam(normalizedcounts, gs)
gsva.es <- gsva(gsvaPar)
gsva.es <- as.data.frame(gsva.es)

gsva.es.t <- t(gsva.es)

gsva.es.table <- as.data.frame(pb.subset@meta.data)
gsva.es.table <- bind_cols(gsva.es.table, gsva.es.t)
gsva.es.table <- gather(gsva.es.table, key = "Geneset", value = "GSVA Score", 53:ncol(gsva.es.table))


#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

gsva.es.table$patientNum <- str_sub(gsva.es.table$patient, start = 4, end = 4)

# Plot GSVA
plotGSVA <- function(metadata, module, t) {
  gData <- filter(gsva.es.table, Geneset == module) 
  gData <- dplyr::filter(gData, time == t)
  gData$treatment <- factor(gData$treatment, levels = c("Placebo", "BAF312", "BAF312-interrupted"))
  lm_fit <- lm(paste(metadata, '~ `GSVA Score`'), gData)
  lm_int <- coef(lm_fit)[1]
  lm_slope <- coef(lm_fit)[2]
  lm_p <- signif(overall_p(lm_fit), digits = 3)
  lm_r <- signif(summary(lm_fit)$r.squared, digits = 3)
  ggplot(gData, aes(x = `GSVA Score`, y = .data[[metadata]], fill = treatment, shape = treatment), color = 'Black') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 12)) +
    labs(title = paste(metadata, "\n", module, "\n", t),
         x = module,
         caption = paste('lm p-val =', lm_p, '; lm r-squared = ', lm_r)) +
    geom_abline(aes(slope = lm_slope, intercept = lm_int), linewidth = 1) +
    scale_fill_viridis(option = "magma", discrete = TRUE) +
    scale_shape_manual(values = 21:23) +
    geom_point(size = 6) +
    geom_text(aes(label = patientNum), size = 5, nudge_x = 0.0, nudge_y = -0.3)
  
  ggsave(paste0("GSVA ", metadata, " ", module, " ", t, " ", myTime(), ".pdf"), height = 4, width = 5)  
}

graphMetadata <- c("Hematoma_Vol_Day1",
                   "PHE_Day1",
                   "Rel_PHE_Day1",
                   "NIHSS_Day1",
                   "MRS_Day90")
graphModules <- unique(gsva.es.table$Geneset)
graphTimes <- c('Day 1', 'Day 3', 'Day 7')

for (met in graphMetadata) {
  for (mod in graphModules) {
    for (t in graphTimes) {
      plotGSVA(met, mod, t)
    }
  }
}









