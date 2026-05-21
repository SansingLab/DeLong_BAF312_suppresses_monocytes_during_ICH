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

mywd <- paste0("[path]/", myTime(), " module scoring")
dir.create(mywd)
setwd(mywd)

timepoints <- c('Day 1', 'Day 3', 'Day 7')


# Load combined downsampled object from script "07 - Combine downsampled objects.R"
load("[path]/s1p.downsampled.RData")

# Functions to plot all timepoints from Placebo and BAF312-treated patients on the same scale
FeaturePlotSingleTimePlacebo_BAF312 <- function(placebo, baf312, feature) {
  minimal <- min(min(placebo[[feature]]), min(baf312[[feature]]))
  maximal <- max(max(placebo[[feature]]), max(baf312[[feature]]))
  ps_placebo <- list()
  ps_baf312 <- list()
  for (t in c('Day 1', 'Day 3', 'Day 7')) {
    placeboSubset <- subset(placebo, subset = time == t)
    p <- FeaturePlot(placeboSubset, features = feature, order = TRUE, pt.size = 0.5) +
      scale_color_viridis_c(limits = c(minimal, maximal), option = "magma") +
      ggtitle(paste("Placebo", t)) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps_placebo[[paste("Placebo", t)]] <- p
  }
  for (t in c('Day 1', 'Day 3', 'Day 7')) {
    baf312Subset <- subset(baf312, subset = time == t)
    p <- FeaturePlot(baf312Subset, features = feature, order = TRUE, pt.size = 0.5) +
      scale_color_viridis_c(limits = c(minimal, maximal), option = "magma") +
      ggtitle(paste("BAF312", t)) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    ps_baf312[[paste("BAF312", t)]] <- p
  }
  ps_list <- list(ps_placebo, ps_baf312)
  return(ps_list)
}

graphModuleScorePlacebo_BAF312 <- function(module) {
  p_list <- FeaturePlotSingleTimePlacebo_BAF312(s1p.downsampled.placebo, s1p.downsampled.baf312, module)
  wrap_plots(p_list[[1]], guides = 'collect', design = "ABC")
  ggsave(paste0("FeaturePlot ordered ModuleScore equal scales for Placebo and BAF312 ", module, " Placebo downsampled ", myTime(), ".pdf"), height = 3, width = 9)
  
  wrap_plots(p_list[[2]], guides = 'collect', design = "ABC")
  ggsave(paste0("FeaturePlot ordered ModuleScore equal scales for Placebo and BAF312 ", module, " BAF312 downsampled ", myTime(), ".pdf"), height = 3, width = 9)
}

# Add module scores of interest
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c("ABCA1", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", 
                                                                     "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", 
                                                                     "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "ACKR3", "CCN1", #"RIGI", 
                                                                     "DENND5A", "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", "EHD1", "EIF1", "ETS2", 
                                                                     "F2RL1", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183", 
                                                                     "HBEGF", "HES1", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5", "IFIH1", "IFIT2", "IFNGR2", "IL12B", "IL15RA", "IL18", "IL1A", 
                                                                     "IL1B", "IL23A", "IL6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", 
                                                                     "KLF9", "KYNU", "LAMB3", "LDLR", "LIF", "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", 
                                                                     "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", "PDE4B", "PDLIM5", 
                                                                     "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLAUR", "PLEK", "PLK2", "PMEPA1", "PNRC1", "PLPP3", "PPP1R15A", "PTGER4", "PTGS2", 
                                                                     "PTPRE", "PTX3", "RCAN1", "REL", "RELA", "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2", "SERPINB8", "SERPINE1", 
                                                                     "SGK1", "SIK1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SPSB1", "SQSTM1", "STAT5A", "TANK", "TAP1", 
                                                                     "TGIF1", "TIPARP", "TLR2", "TNC", "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9", "TNFSF9", "TNIP1", "TNIP2", "TRAF1", 
                                                                     "TRIB1", "TRIP10", "TSC22D1", "TUBB2A", "VEGFA", "YRDC", "ZBTB10", "ZC3H12A", "ZFP36")), 
                                  name = "TNFA_SIGNALING_VIA_NFKB", search = TRUE)

# Genes upregulated at least 2-fold in Bai 2024, ICH vs Hypertensive
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c("NECAB1", "FABP2", "CD177", "ADAMTS2", "OLAH", "DCC", "SLC1A3", "EPB41L4B", "MMP8", 
                                                                     "RP11-752G15.4", "METTL7B", "AMPH", "GPR84", "HP", "GPR64", "TGM5", "C19orf59", 
                                                                     "S100A12", "DEFA8P", "PCOLCE2", "SLC51A", "AC006483.5", "C1orf192", "PRL", 
                                                                     "TNFAIP8L3", "ANKRD22", "EDNRB", "NEBL", "ARG1", "AC004069.1", "RNA5SP180", "HPD", 
                                                                     "KB-1517D11.2", "ANXA3", "SNORD74", "IGLC7", "RP11-497H16.4", "TSKU", "AC010731.6", 
                                                                     "AC092299.6", "PTGFR", "PPARG", "MMP9", "SLC9A4", "C6orf195", "CRISP3", "KAL1", 
                                                                     "CLEC5A", "FAM20A", "CCDC140", "DPRXP2", "AP3B2", "HIST1H2AM", "GVINP2", "CYP19A1", 
                                                                     "CEACAM8", "AMPD1", "OLR1", "RNASE1", "BMX", "ITGA8", "OVCH1", "TDRD9", "IL10", 
                                                                     "LRRN1", "RP11-443P15.2", "S100A8", "OLFM4", "NSG2", "CTB-50E14.5", "ADAMTS3", 
                                                                     "FCGR1C", "PROK2", "SDC1", "SERPINB10", "NDNF", "TMEM52B", "SLCO2B1", "IGHA2", 
                                                                     "CEACAM6", "FCGR1A", "AC005795.1", "KCNMA1", "DEFA3", "RP11-497H16.6", "AC073342.1", 
                                                                     "RP11-925D8.1", "FAM124A", "AL353791.1", "KL", "ZDHHC19", "ZNF608", "BCL2A1", 
                                                                     "COL17A1", "PRSS12", "CLEC4D", "ATP2C2", "RP11-277J24.1", "IGHG2", "ORM1", 
                                                                     "HTRA1", "OOSP1", "MMP2", "TLR5", "INHBA", "GLDN", "CEACAM1", "C15orf65", 
                                                                     "RP11-925D8.3", "RPEP6", "SCN9A", "IGLC3", "IGJ", "C6orf223", "AC020636.1", 
                                                                     "RNASE2", "PROKR2", "SMIM10", "RP11-1220K2.2", "OGFOD1P1", "TCN1", "LTF", 
                                                                     "ANKRD34B", "KIAA1244", "NDST3", "RN7SKP26", "DDX10P1", "RP11-111G23.1", "TMTC1", 
                                                                     "RP11-925D8.6", "C11orf82", "SFRP1", "SLC7A11", "PTGR1", "OMG", "FCGR1B", "ERG", 
                                                                     "ABCA13", "SETP11", "TJP1", "RP11-497H16.2", "NACAP3", "STAB2", "RP11-1319K7.1", 
                                                                     "RP11-925D8.5", "GYG1", "WDFY3-AS2", "PFKFB2", "ASPH", "F5", "CAV1", "SLC26A8", 
                                                                     "CA12", "VSIG4", "KLF5", "IRAK3", "RSL24D1P6", "Y_RNA", "snoU13", "TCTEX1D1", 
                                                                     "RN7SKP97", "CAMP", "FIGN", "AC132872.2", "ATP9A", "SEPT7P8", "IGLV7-43", 
                                                                     "FAM106DP", "NAIP", "SMPDL3A", "CRISP2", "SHOX2", "NAMPTL", "RNU6-37P", 
                                                                     "SPATA20P1", "C7", "MIXL1", "B3GNT5", "CASP5", "PDE6H", "IL1R2", "SYN2", "SLPI", 
                                                                     "RP1-34B20.4", "CCDC36", "GRB10", "AC112200.1", "SULT1B1", "SAMD15", "CKAP4", 
                                                                     "SEMA6A", "UBE2J1", "UCHL1", "RP11-1023L17.2", "IGKV2-30", "MANSC1", "VNN1", 
                                                                     "RNU7-29P", "RP11-98J23.1", "RNU6-1197P", "C8orf88", "HPGD", "ITGA7", "PLBD1", 
                                                                     "MGST1", "QPCT", "TMEM45A", "PYGL", "PSTPIP2", "CNTNAP3B", "LARP1P1", "MYO10", 
                                                                     "MAOB", "PLAU", "PGLYRP1", "HGF", "IGHA1", "ALDH1L2", "IGHV6-1", "SCRG1", "LCN2", 
                                                                     "LSMEM1", "SERPINB2", "CCDC73", "RPS3AP43", "RPL30P7", "RPH3A", "S100A9", 
                                                                     "DACH1", "STAC", "SLC37A3", "GADD45A", "BPI", "CAPN13", "RP11-632F7.1", "INSC", 
                                                                     "TLR8", "ZNF438", "HIST1H2BM", "RP11-442P12.2", "SAMSN1", "NOS1AP", "IGLV1-40", 
                                                                     "CARD17", "UGCG", "CES1", "CYSTM1", "HMMR", "IL18RAP", "TUBBP5", "RGL4", 
                                                                     "RPS3AP32", "STOX2", "C1orf226", "BCL6", "PNPLA1", "AC138123.2", "NRN1", "IL18R1", 
                                                                     "FGF13", "LMNB1", "RN7SL271P", "KLHL2", "PLSCR1", "MS4A4A", "GSDMC", "FOXC1", 
                                                                     "RNU6-176P", "RP11-848G14.5", "NFIL3", "Metazoa_SRP", "RP11-420K14.1", "AGPAT9", 
                                                                     "AC000081.2", "S100P", "FOLR3", "CTD-2330J20.1", "RNU6-80P", "CPE", "RPL7P24", 
                                                                     "LY6G6C", "TMEM253", "SIPA1L2", "SLC2A3", "PRDM5", "ENTPD7", "GNG10", "CARD6", 
                                                                     "FCER1G", "OR9A2", "HSPA4L", "PPIAP29", "EFCAB2", "MS4A3", "SRPK1", "MAP2K6", 
                                                                     "FKBP5", "OR2B6", "DLC1", "ADAM9", "C5orf47", "RP11-95K23.7", "FCAR", "SELL", 
                                                                     "TPST1", "A2M", "MAB21L3", "RCVRN", "RNU6-1013P", "BEND7", "HNRNPA1P52", 
                                                                     "HIST2H2BE", "IGLC2", "IGKV3-20", "MGAM", "TIMP3", "PDLIM1P4", "NKAIN2", "SUCNR1", 
                                                                     "FAM151B", "LYVE1", "GLDC", "CDC6", "HORMAD1", "FGD4", "SIGLECL1", 
                                                                     "RP11-25K21.6", "TEAD3", "HTATSF1P2", "CEP55", "CR1", "C5orf27", "CYP1B1", 
                                                                     "RNU6-1300P", "SPINK8", "RNU6-313P", "MED6P1", "ZNF215", "RP11-65E22.2", 
                                                                     "ALOX5AP", "RAI2", "DUSP13", "MZB1", "CNTNAP3", "MLTK", "MND1", "AC023050.1", 
                                                                     "NCAPG", "AQP9", "MAPK14", "OR7E140P", "AC022431.2", "MIR5690", "OR6N2", "CEP19", 
                                                                     "RNU6-574P", "FLT3", "GAPDHP33", "ECHDC3", "TRIM6", "LIPN", "NEK2", "PGD", 
                                                                     "SPCS2P4", "NLRC4", "FAM169B", "SLC8A1", "XKR3", "BASP1", "WDFY3", "CKAP2L")),
                                  name = 'Bai_2024_ICHvHTN_Up', search = TRUE)

# Genes downregulated at least 2-fold in Bai 2024, ICH vs Hypertensive
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c("TRAJ44", "SLC22A18", "TMEM204", "ZFP92", "TTC16", "RP11-108K14.4", "SIGLEC17P", 
                                                                     "CHST8", "COL13A1", "MIR342", "ACRBP", "RP11-1055B8.7", "LTK", "KLRK1", "MEIS3P2", 
                                                                     "FCRL6", "MIR4446", "GGN", "TTYH3", "WNT10A", "KRT2", "KCNK10", "ZAP70", "IL9R", 
                                                                     "PDGFRB", "EPHA2", "RN7SL5P", "GAA", "HIST1H3H", "ZNF683", "RP11-553L6.3", 
                                                                     "VIPR2", "GLB1L2", "GZMH", "GAS2L1", "R3HDM4", "LGR6", "S1PR5", "RN7SL838P", 
                                                                     "MIR1276", "ERBB2", "KCNG1", "AC004017.1", "GNLY", "SEMA4C", "TSEN54", 
                                                                     "RP13-39P12.2", "ADAMTS10", "CDR2L", "PDZD4", "DUSP2", "RN7SL653P", "CALHM1", 
                                                                     "MIR661", "OLA1P2", "NPTX1", "RN7SL859P", "PODN", "CACNA1H", "DEGS2", "SCN11A", 
                                                                     "RN7SL67P", "CACNA2D2", "AC012314.19", "WNT10B", "AC079781.8", "MIR3194", 
                                                                     "COL6A2", "ENHO", "KCNT1", "IGLV3-21", "NOG", "CCL4L1", "IFI27", "EMILIN1", 
                                                                     "LLGL2", "AC022400.2", "MIR3142", "MYRF", "CARNS1", "MIR4637", "PTCHD2", "CRYBB3", 
                                                                     "LGALS9B", "CSDC2", "KCNA5", "PTCRA", "RP11-436I9.5", "PTPRN", "PTGDS", "LYL1", 
                                                                     "EBF4", "SLC1A7", "VSTM2B", "MIR4473", "BZRAP1", "MIR1250", "TRBJ2-2", "GPR56", 
                                                                     "CCDC79", "LIM2", "MIR4772", "SEPT1", "CABP4", "WNT1", "SMPD4P1", "SPON2", "BAI2", 
                                                                     "MIRLET7E", "RNF208", "PPT2-EGFL8", "LOC440461", "MATK", "AQP7", "MXRA8", "KIF19", 
                                                                     "MIR181A2", "GP9", "CLIC3", "ALG1L11P", "HBQ1", "PRSS30P", "NYAP1", "GSC", "SBK1", 
                                                                     "RP11-229P13.19", "C8G", "AC093788.1", "AC092566.1", "NOXO1", "RN7SL842P", 
                                                                     "AC108456.1", "ABHD17AP6", "TRBJ2-7", "LCNL1", "CLDND2", "MIR150", "GZMM", 
                                                                     "MIR4669", "RP11-366O20.5", "SHROOM2", "HIST1H3F", "MIR3687", "SCT", "MIR4439", 
                                                                     "MIR33A", "MIR3648")),
                                  name = 'Bai_2024_ICHvHTN_Down', search = TRUE)

# Of the top 100 DEGs in Walsh et al swine, these are upregulated
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c('FN1', 'GJC1', 'MAPK10', 'HTR7', 'NETO2', 'ADAMTS2', 'NT5DC2', 'P2RY2',
                                                                     'CNGA4', 'CRISPLD2', 'COL4A4', 'C4BPA', 'ADORA1', 'MYBPH', 'DGAT2', 'CD14',
                                                                     'MEGF9', 'MCEMP1', 'TINAGL1', 'TBC1D8B', 'NRG1', 'TGM1', 'TFPI2', 'TGM3',
                                                                     'IL1R2', 'CHIT1', 'S100A12', 'S100A9', 'S100A8', 'LINGO3', 'LTF', 'RETN',
                                                                     'OLIG1', 'TNFAIP6', 'TMEM150C', 'IL34', 'FOS', 'CXCL8', 'CXCL2', 'EPB41L3',
                                                                     'EREG', 'IL18', 'PTGS2', 'SEMA3C', 'VSIG4', 'CD163', 'ADGRE1', 'KLHL23',
                                                                     'TCEA3', 'WFDC1', 'SFRP5', 'TMEM26', 'NME3', 'MS4A7')),
                                  name = 'Walsh_2019_Up', search = TRUE)

# Of the top 100 DEGs in Walsh et al (swine), these are downregulated
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c('SERPINE3', 'DCHS1', 'FNDC4', 'CNIH3', 'L1CAM', 'NPSR1', 'TMEM196', 'SLC45A2',
                                                                     'IL33', 'BOC', 'DOC2A', 'CASKIN1', 'SPDEF', 'CEACAM19', 'BPIFB4', 'ZKSCAN2',
                                                                     'LEXM', 'MGAT5B', 'FAAH', 'GPR179', 'USP51', 'NTNG2', 'ATP4A', 'IFN-DELTA-3',
                                                                     'FAM171A1', 'CLEC2L', 'SCIN', 'MICU3', 'IFN-DELTA-7', 'HPD', 'KRT2', 'ANK2',
                                                                     'DYSF', 'WEE2', 'TP53I3', 'CACNG1', 'ZNF835', 'MYH2', 'IFN-DELTA-9', 'HTR4',
                                                                     'SEMA5B', 'CSMD2', 'TMEM176A', 'IFN-ALPHA-4', 'C16orf90', 'RABGGTA')),
                                  name = 'Walsh_2019_Down', search = TRUE)

# Genes upregulated at least 1.6-fold in ICH vs healthy control - Durocher 2019
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c("GYG2-AS1", "MAPK14", "RP3-523E19.2", "RNU6-1254P", "GK3P", "ACSL4", "IRAK3", "MCTP2", "LINC00211", "IL1RAP", 
                                                                     "PHTF1", "GNG10", "RNU6-1052P", "GK-AS1", "MGAM", "SH3GLB1", "WDFY3", "ST3GAL4", "FKBP5", "RNU6-96P", "CR1", 
                                                                     "ZNF438", "RNU6-1199P", "RP11-242C19.2", "AP001434.2", "RNU6-1054P", "SULT1B1", "RP1-55C23.7", "RNU6-241P", 
                                                                     "RP11-597D13.9", "RNU6-1019P", "IL4R", "NLRC4", "ASPH", "MIR618", "LRG1", "CASP5", "IL1R2", "NAIP", "CEBPB", 
                                                                     "FCAR", "SIPA1L2", "TLR5", "RNU6-1003P", "ACSL1", "DYSF", "LMNB1", "LINC00694", "F5", "CD59", "CYP1B1", "NFIL3", 
                                                                     "RNU6-1152P", "AC007278.3", "DSC2", "RP11-585F1.2", "RNU7-163P", "MIR21", "RP11-1198D22.2", "SLC37A3", "UGCG", 
                                                                     "SLC26A8", "CD163", "KREMEN1", "RP11-701P16.1", "RNU6-705P", "PLSCR1", "FCGR1A", "KLHL2", "PROK2", 
                                                                     "XX-C2158C6.1", "FCGR1C", "BCL2A1", "ALPL", "RP11-76E17.3", "IL18R1", "RNU6-226P", "CYSTM1", "GYG1", "PFKFB3", 
                                                                     "ORM2", "IL18RAP", "BMX", "ANXA3", "GYG1P1", "CLEC4D", "CD177P1", "C19orf59",
                                                                     "RNF24", "EXOC6", "RNU7-16P", "FAR1P1", "RP11-537H15.4", "KCNJ2", "LINC01093", "GNAQP1", "RGS2", "NRBF2P2", 
                                                                     "AL139812.1", "PPP4R1", "RNU7-126P", "EGLN1P1", "AL354933.1", "RP5-968J1.1", "MSL1", "RNU6-674P", "APMAP", 
                                                                     "ATP6V1C1", "IDI1", "LILRA6", "ERLIN1", "AGFG1", "SLC2A14", "HLX", "FOSL2", "DNAJC3", "SLC2A3P4", "LINC01094", 
                                                                     "ATP11B", "MIR4476", "RNU6-797P", "RP11-1334A24.6", "STK3", "DHRS13", "MSRB1P1", "RNU7-11P", "RP11-585F1.6", 
                                                                     "TMPRSS7", "GLT1D1", "PIK3AP1", "MIR617", "MMP25", "RP11-217B1.2", "RBMS1P1", "VCAN", "SLC2A3P2", 
                                                                     "RP11-326A19.5", "TP53I11", "SIGLEC9", "LTBR", "IMPDH1P4", "RP11-99E15.2", "TLR2", "IL1R1", "AC007743.1", 
                                                                     "OR8G2P", "RNU7-9P", "IFNAR1", "LINC00265-3P", "AL445665.2", "MIR4802", "WSB1", "CPD", "RNU7-47P", "TLR4", 
                                                                     "CLIC1P1", "GADD45A", "KIF1B", "FUT7", "LINC00266-3", "CHSY1", "CFLAR-AS1", "GCA", "RNU6-917P", 
                                                                     "RP1-232L22__B.1", "ADAM9", "TRIM25", "LINC00999", "RNU6-1308P", "FAM157A", "RNU6-587P", "DGAT2", "TBC1D8", 
                                                                     "FGD4", "RP11-295G20.2", "ETS2", "RP11-274B21.4", "MAP2K6", "GRAMD1A", "RNU6-725P", "PGS1", "RNU6-344P", 
                                                                     "RNU6-855P", "TRIQK", "RP11-473O4.4", "RNU6-430P", "SERPINB1", "HK3", "MIR3614", "CR1L", "CNEP1R1", 
                                                                     "RNU6-1146P", "SLC2A3", "RP11-504P24.2", "LIMK2", "CDC42EP3", "GALNT14", "LIN7A", "SAMSN1", "CEACAM1", 
                                                                     "FAM41C", "MRVI1", "AC002511.2", "MERTK", "RNU6-567P", "MEF2AP1", "B4GALT5", "MKNK2P1", "FFAR2", 
                                                                     "RP11-153M7.5", "PYGL", "CA4", "RP1-177A13.1", "MAK", "UBE2J1", "FAM126B", "B3GNT5", "RNU6-719P", 
                                                                     "CTB-167B5.2", "VENTXP2", "MPZL3", "HRH2", "RP11-291L22.4", "ROPN1L", "ZNF608", "LINC00266-4P", "MIR4441", 
                                                                     "GAS7", "EXT1", "RP11-1280N14.3", "MKNK1", "RNU6-177P", "GK", "FLOT2", "GPR97", "LILRA5", "GPR141", "SRPK1", 
                                                                     "RNU6-196P"
                                                                     )),
                                  name = 'Durocher_2019_Up', search = TRUE)

# Genes downregulated at least 1.6-fold in ICH vs healthy control - Durocher 2019
s1p.downsampled <- AddModuleScore(s1p.downsampled, features = list(c("TRBV20-1", "TRAJ26", "TRAJ45", "TRAJ22", "TRAJ13", "KLRB1", "TRAJ27", "TRAJ32", "TRAJ41", "TRAJ42", "TSTA3", 
                                                                     "TRAJ35", "CD96", "TRAJ36", "TRAJ9", "TRBV20OR9-2", "AC087499.4", "TRAJ10", "TRAV9-2", "CD2", "TRAJ38", "GBP4", 
                                                                     "TESPA1", "VSIG1", "TRAV4", "TRAJ31", "TRAC", "TRAJ37", "TRAJ1", "TRAV12-3", "TRBV3-1", "TRAJ16", "THEMIS", 
                                                                     "CD27", "GPR174", "CTC-505O3.3", "TRAJ14", "TRAJ34", "CD3G", "TRAV16", "STAT4", "SNORD108", "TRAJ44", "IL7R", 
                                                                     "TRAJ7", "C14orf1", "TRAJ30", "TRAJ3", "TRAJ18", "TRBC2", "ANKRD36P1", "TRAJ48", "CD40LG", "MAL", "RP11-58E21.5", 
                                                                     "TRAJ29", "TRAJ40", "TRAJ39", "TRAJ20", "SLAMF6", "TRAJ2", "TRAJ11", "CD28", "TRAJ19", "LDHB", "CTD-2049O17.1", 
                                                                     "RASGRP1", "SLFN12L", "ITK", "PRKCH", "PCED1B-AS1", "SLC38A5", "OPTN", "RP11-432F4.2", "TC2N", "IKZF3", 
                                                                     "RP11-1036F1.1", "DPP4", "NFATC2", "CTC-523E23.11", "CBLB", "LINC00861", "TRAJ6", "CCND2", "RP11-63K6.4", 
                                                                     "PRKCQ-AS1", "TRAJ21", "OR7E38P", "CCL5", "PARP15", "NELL2", "RP11-697N18.1", "CD3E", "RPLP2P1", "RORA", 
                                                                     "RP11-554I8.1", "RP3-324O17.4", "TRAJ50", "RP1-140K8.2", "AC107983.4", "TRAJ17", "RP11-301G21.1", "C14orf64", 
                                                                     "TRBV7-6", "INPP4B", "TRBV6-5", "RPL37P15", "TRAV8-3", "TRBJ2-1", "ARL4C", "NPM1P6", "RP11-155G14.1", 
                                                                     "SLC38A1", "RP1-97D16.1", "C9orf43", "SIDT1", "ISCUP1", "TRAJ23", "UBASH3A", "ABLIM1", "SEPT1", "OXNAD1", 
                                                                     "RASGRF2", "LCK", "CTB-47B8.1", "DOCK9", "SLAMF1", "PRKCQ", "RP11-1072N2.2")),
                                  name = 'Durocher_2019_Down', search = TRUE)


s1p.downsampled.placebo <- subset(s1p.downsampled, treatment == "Placebo")
s1p.downsampled.baf312 <- subset(s1p.downsampled, treatment == "BAF312")

save.image("module_scoring.RData")

load('[path]/module_scoring.RData')


# UMAP graphs of TNFa signaling
graphModuleScorePlacebo_BAF312("TNFA_SIGNALING_VIA_NFKB1")

# Summary dotplots of DEGs from previous studies
modules <- c('Bai_2024_ICHvHTN_Up1',
             'Durocher_2019_Up1', 
             'Walsh_2019_Up1',
             'Bai_2024_ICHvHTN_Down1',
             'Durocher_2019_Down1',
             'Walsh_2019_Down1'
             )
modules <- factor(modules, levels = modules)

clusters <- c('Monocytes', 'CD16+ Monocytes', 'NK cells', 'CD4+ T cells', 'CD8+ T cells', 'B cells')

module.means <- data.frame(Cluster = character(),
                           Timepoint = character(),
                           Module = character(),
                           Direction = character(),
                           `Module Score` = numeric())

for (c in clusters) {
  s1p.c <- subset(s1p.downsampled.placebo, lineageClusterIdents == c)
  for (t in timepoints) {
    s1p.t <- subset(s1p.c, time == t)
    for (m in modules) {
      temp <- data.frame(Cluster = c,
                         Timepoint = t,
                         Module = m,
                         Direction = ifelse(grepl('Down', m), 'Down', 'Up'),
                         `Module Score` = mean(s1p.t@meta.data[[m]])
                         )
      module.means <- bind_rows(module.means, temp)
    }
  }
}

module.means$Module <- factor(module.means$Module, levels = modules)
module.means$Direction <- factor(module.means$Direction, levels = c('Up', 'Down'))
module.means$Cluster <- factor(module.means$Cluster, levels = rev(clusters))


# Facet by direction of regulation
ggplot(module.means, aes(x = Module, y = Cluster, color = Module.Score)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(size = 8) +
  scale_color_viridis_c(option = "magma") +
  facet_wrap(.~Direction, scales = 'free_x')
ggsave(paste0('Module Score ICH vs healthy PB summary merged timepoints facetted by Direction ', myTime(), '.pdf'), height = 4, width = 5.5)

# Facet by direction of regulation - Day 1
filter(module.means, Timepoint == "Day 1") %>%
  ggplot(aes(x = Module, y = Cluster, color = Module.Score)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_point(size = 8) +
    scale_color_viridis_c(option = "magma") +
    facet_wrap(.~Direction, scales = 'free_x')
ggsave(paste0('Module Score ICH vs healthy PB summary day 1 facetted by Direction ', myTime(), '.pdf'), height = 4, width = 5.5)


