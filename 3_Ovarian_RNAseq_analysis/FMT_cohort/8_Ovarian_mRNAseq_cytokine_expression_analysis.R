options(stringsAsFactors = FALSE)

# load libraries for analysis
library(pheatmap)
library(matrixStats)
library(beeswarm)
library(outliers)

##############################################
# Menopause-microbiome project
# Generate heatmap of cytokines
##############################################

##############################################
# 1. Load count data
##############################################

my.counts <- read.table("./2024-03-29_MM_FMT_Ovaries_DESeq2_model_with_FMT__log2_counts_matrix_DEseq2_SVA.txt", header = TRUE)

##############################################
# 2. Generate heatmap
##############################################

cytokine_genes <- c(
  "Tnf", "Il6", "Il12a", "Il12b", "Il14","Il17f", "Ccl2", "Cxcl2", "Cxcl9",
  "Cxcl10", "Nlrp3", "Il10", "Il4", "Il13", "Tgfb1", "Tgfb2"
)

present_genes <- intersect(cytokine_genes, rownames(my.counts))
data.for.heatmap <- my.counts[present_genes, ]

row_var <- rowVars(as.matrix(data.for.heatmap))
data.filtered <- data.for.heatmap[row_var > 0, ]

pdf(paste0(Sys.Date(), "_FMT_mRNAseq_cytokine_heatmap.pdf"))
pheatmap(data.filtered,
         cluster_cols = F,
         cluster_rows = F,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T, scale="row",
         main = "FMT mRNAseq - cytokines", cellwidth = 15)
dev.off()

################################################################################
sink(file = paste(Sys.Date(),"MM_Ovarian_aging_FMT_cytokine_expression_analysis_session_Info.txt",sep="_"))
sessionInfo()
sink()
