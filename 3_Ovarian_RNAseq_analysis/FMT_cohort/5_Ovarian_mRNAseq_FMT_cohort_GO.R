setwd('/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Github/Ovarian_mRNAseq/FMT_cohort/GSEA/')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(phenoTest)
library(qusage)
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

library(clusterProfiler)
library(enrichplot)

set.seed(123123)

###############################################
# Menopause-Microbiome project
# Ovarian mRNA-seq - FMT cohort (FMT-YF vs. FMT-EF)
# Run GSEA
###############################################

# Load GSEA R functions and datasets
source("./2_Ovarian_mRNAseq_FMT_cohort_GSEA_function.R")

load("./GSEA_DBs/2024-03-29_GeneSetCollections_for_Phenotest_GSEA_mouse.RData")
load("./GSEA_DBs/2024-03-29_MM_FMT_Expressed_TF_summary_GeneSetCollections_for_Phenotest_GSEA.RData")

# Load DESeq2 data
load("../DESeq2/2024-03-29_MM_FMT_Ovaries_FMT.RData")

##### Prepare gene list using DESeq2 t-statistic to rank genes #####

my.MM.FMT.ovary <- data.frame(res.FMT)
my.MM.FMT.ovary$Gene.name <- rownames(my.MM.FMT.ovary)

MM.FMT.ovary.genelist <- my.MM.FMT.ovary$stat
names(MM.FMT.ovary.genelist) <- my.MM.FMT.ovary$Gene.name
MM.FMT.ovary.genelist <- sort(MM.FMT.ovary.genelist, decreasing = TRUE)

##### Run GSEA #####

MM.FMT.ovary.ENS.GO.BP      <- run_enrich(MM.FMT.ovary.genelist, "MM_Ovary_FMT_effect", my.fdr = 0.05, my.ontology = Sym.ENS.GO.BP, my.ontology.name = "ENS.GO.BP")
MM.FMT.ovary.ENS.GO.CC      <- run_enrich(MM.FMT.ovary.genelist, "MM_Ovary_FMT_effect", my.fdr = 0.05, my.ontology = Sym.ENS.GO.CC, my.ontology.name = "ENS.GO.CC")
MM.FMT.ovary.ENS.GO.MF      <- run_enrich(MM.FMT.ovary.genelist, "MM_Ovary_FMT_effect", my.fdr = 0.05, my.ontology = Sym.ENS.GO.MF, my.ontology.name = "ENS.GO.MF")

MM.FMT.ovary.TF.targets     <- run_enrich(MM.FMT.ovary.genelist, "MM_Ovary_FMT_effect", my.fdr = 0.1, my.ontology = Sym.m3.gtrd, my.ontology.name = "TF_GTRD")

##### Plot Bubble plots #####

plot_bubble_plot("MM_FMT_ovary_ENS.GO.BP", "./2024-03-29_MM_Ovary_FMT_effect_ENS.GO.BP_FDR_5_Phenotest_GSEA_Analysis_table_485_significant.txt")
plot_bubble_plot("MM_FMT_ovary_ENS.GO.CC", "./2024-03-29_MM_Ovary_FMT_effect_ENS.GO.CC_FDR_5_Phenotest_GSEA_Analysis_table_80_significant.txt")
plot_bubble_plot("MM_FMT_ovary_ENS.GO.MF", "./2024-03-29_MM_Ovary_FMT_effect_ENS.GO.MF_FDR_5_Phenotest_GSEA_Analysis_table_102_significant.txt")

plot_bubble_plot("MM_FMT_ovary_TF_GTRD", "./2024-03-29_MM_Ovary_FMT_effect_TF_GTRD_FDR_10_Phenotest_GSEA_Analysis_table_74_significant.txt")


# SenMayo

my.SEN.MAYO <- Sym.SEN.MAYO$SAUL_SEN_MAYO

TERM2GENE <- data.frame(TERM = rep("SENMAYO", length(my.SEN.MAYO)), 
                        GENE = my.SEN.MAYO, 
                        stringsAsFactors=FALSE)

gseaResult.SENMAYO <- GSEA(MM.FMT.ovary.genelist, 
                           TERM2GENE = TERM2GENE, 
                           minGSSize = 15, 
                           maxGSSize = 500, 
                           eps = 1e-100,
                           pvalueCutoff = 1,
                           pAdjustMethod = 'BH',
                           by = 'fgsea',
                           seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_SEN_MAYO.pdf"))
gseaplot2(gseaResult.SENMAYO, geneSetID = c("SENMAYO"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()

# Up with age - "Induces"

TERM2GENE <- data.frame(TERM = rep("CellAge_UP", length(my.cellage.upwithage$'Gene symbol')), 
                        GENE = my.cellage.upwithage$'Gene symbol', 
                        stringsAsFactors=FALSE)

gseaResult.CELLAGE.UP <- GSEA(MM.FMT.ovary.genelist, 
                           TERM2GENE = TERM2GENE, 
                           minGSSize = 15, 
                           maxGSSize = 500, 
                           eps = 1e-100,
                           pvalueCutoff = 1,
                           pAdjustMethod = 'BH',
                           by = 'fgsea',
                           seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_Cell_Age_UPWITHAGE.pdf"))
gseaplot2(gseaResult.CELLAGE.UP, geneSetID = c("CellAge_UP"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()


# Down with age - "Inhibits"

TERM2GENE <- data.frame(TERM = rep("CellAge_DOWN", length(my.cellage.downwithage$'Gene symbol')), 
                        GENE = my.cellage.downwithage$'Gene symbol', 
                        stringsAsFactors=FALSE)

gseaResult.CELLAGE.DOWN <- GSEA(MM.FMT.ovary.genelist, 
                              TERM2GENE = TERM2GENE, 
                              minGSSize = 15, 
                              maxGSSize = 500, 
                              eps = 1e-100,
                              pvalueCutoff = 1,
                              pAdjustMethod = 'BH',
                              by = 'fgsea',
                              seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_Cell_Age_UPWITHAGE.pdf"))
gseaplot2(gseaResult.CELLAGE.DOWN, geneSetID = c("CellAge_DOWN"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()


# Plot expression levels of Ncoa1 and Usp7

load("../DESeq2/2024-03-29_MM_FMT_Ovaries_DESeq2_Analysis_STAR_post_SVA_normalized_counts.RData")

my.colors <- rep("#9C509F",14)
my.colors[grep("FMT_EF",colnames(tissue.cts))] <- "#F99E26"

pdf(paste(Sys.Date(),"_Normalized_counts_Ncoa1_expression_barplot.pdf"))
barplot(tissue.cts["Ncoa1",], ylab = "Normalized log2(counts) Ncoa1 expression", las = 2, col = my.colors)
dev.off()

pdf(paste(Sys.Date(),"_Normalized_counts_Usp7_expression_barplot.pdf"))
barplot(tissue.cts["Usp7",], ylab = "Normalized log2(counts) Usp7 expression", las = 2, col = my.colors)
dev.off()

################################################################################
sink(file = paste(Sys.Date(), "MM_FMT_GSEA_session_Info.txt", sep ="_"))
sessionInfo()
sink()
