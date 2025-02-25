setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Github/Ovarian_mRNAseq/FMT_cohort/GSEA")
options(stringsAsFactors = F)

library(clusterProfiler)
library(DESeq2)
library(enrichplot)

###############################
# Menopause-microbiome project - FMT RNAseq analysis
# Candidate TFs GSEA analysis
# Peak info downloaded from cistrome.org
###############################

###############################
# 1. Import data
###############################

# Load DESeq2 data
load("../DESeq2/2024-03-29_MM_FMT_Ovaries_FMT.RData")

##### Prepare gene list using DESeq2 t-statistic to rank genes #####

my.MM.FMT.ovary <- data.frame(res.FMT)
my.MM.FMT.ovary$Gene.name <- rownames(my.MM.FMT.ovary)

MM.FMT.ovary.genelist <- my.MM.FMT.ovary$stat
names(MM.FMT.ovary.genelist) <- my.MM.FMT.ovary$Gene.name
MM.FMT.ovary.genelist <- sort(MM.FMT.ovary.genelist, decreasing = TRUE)

# Annotated TF target gene lists

my.Foxo1.genes  <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Foxo1_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Nucks1.genes <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Nucks1_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Ikzf1.genes  <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Ikzf1_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Ing1.genes   <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Ing1_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Lyl1.genes   <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Lyl1_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Mllt3.genes  <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Mllt3_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Ncoa1.genes  <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Ncoa1_ChIPseq_HOMER_tss_peaks_gene_names.txt")
my.Usp7.genes   <- read.table("../TF_peaks_files/2023-10-12_Downloaded_Usp7_ChIPseq_HOMER_tss_peaks_gene_names.txt")

###############################
# 2. Run GSEA
###############################

set.seed(123123)

# Foxo1

TERM2GENE <- data.frame(TERM = rep("Foxo1", length(my.Foxo1.genes)), 
                        GENE = my.Foxo1.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Foxo1 <- GSEA(MM.FMT.ovary.genelist, 
                         TERM2GENE = TERM2GENE, 
                         minGSSize = 15, 
                         maxGSSize = 500, 
                         eps = 1e-100,
                         pvalueCutoff = 1,
                         pAdjustMethod = 'BH',
                         by = 'fgsea',
                         seed = TRUE)

# no term enriched under specific pvalueCutoff...

# Nucks1

TERM2GENE <- data.frame(TERM = rep("Nucks1", length(my.Nucks1.genes)), 
                        GENE = my.Nucks1.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Nucks1 <- GSEA(MM.FMT.ovary.genelist, 
                          TERM2GENE = TERM2GENE, 
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          eps = 1e-100,
                          pvalueCutoff = 1,
                          pAdjustMethod = 'BH',
                          by = 'fgsea',
                          seed = TRUE)

# no term enriched under specific pvalueCutoff...

# Nucks1 published peaks

TERM2GENE <- data.frame(TERM = rep("Nucks1_p", length(my.Nucks1.peak.genes.published$Gene)), 
                        GENE = my.Nucks1.peak.genes.published$Gene, 
                        stringsAsFactors=FALSE)

gseaResult.Nucks1p <- GSEA(MM.FMT.ovary.genelist, 
                           TERM2GENE = TERM2GENE, 
                           minGSSize = 15, 
                           maxGSSize = 500, 
                           eps = 1e-100,
                           pvalueCutoff = 1,
                           pAdjustMethod = 'BH',
                           by = 'fgsea',
                           seed = TRUE)

# no term enriched under specific pvalueCutoff...

# Ikzf1

TERM2GENE <- data.frame(TERM = rep("Ikzf1", length(my.Ikzf1.genes)), 
                        GENE = my.Ikzf1.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Ikzf1 <- GSEA(MM.FMT.ovary.genelist, 
                         TERM2GENE = TERM2GENE, 
                         minGSSize = 15, 
                         maxGSSize = 500, 
                         eps = 1e-100,
                         pvalueCutoff = 1,
                         pAdjustMethod = 'BH',
                         by = 'fgsea',
                         seed = TRUE)

# no term enriched under specific pvalueCutoff...

# Ing1

TERM2GENE <- data.frame(TERM = rep("Ing1", length(my.Ing1.genes)), 
                        GENE = my.Ing1.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Ing1 <- GSEA(MM.FMT.ovary.genelist, 
                        TERM2GENE = TERM2GENE, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        eps = 1e-100,
                        pvalueCutoff = 1,
                        pAdjustMethod = 'BH',
                        by = 'fgsea',
                        seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_Ing1_FDR100.pdf"))
gseaplot2(gseaResult.Ing1, geneSetID = c("Ing1"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()

# Lyl1

TERM2GENE <- data.frame(TERM = rep("Lyl1", length(my.Lyl1.genes)), 
                        GENE = my.Lyl1.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Lyl1 <- GSEA(MM.FMT.ovary.genelist, 
                        TERM2GENE = TERM2GENE, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        eps = 1e-100,
                        pvalueCutoff = 1,
                        pAdjustMethod = 'BH',
                        by = 'fgsea',
                        seed = TRUE)

# no term enriched under specific pvalueCutoff...

# Mllt3

TERM2GENE <- data.frame(TERM = rep("Mllt3", length(my.Mllt3.genes)), 
                        GENE = my.Mllt3.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Mllt3 <- GSEA(MM.FMT.ovary.genelist, 
                         TERM2GENE = TERM2GENE, 
                         minGSSize = 15, 
                         maxGSSize = 500, 
                         eps = 1e-100,
                         pvalueCutoff = 1,
                         pAdjustMethod = 'BH',
                         by = 'fgsea',
                         seed = TRUE)

# no term enriched under specific pvalueCutoff...

# Ncoa1

TERM2GENE <- data.frame(TERM = rep("Ncoa1", length(my.Ncoa1.genes)), 
                        GENE = my.Ncoa1.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Ncoa1 <- GSEA(MM.FMT.ovary.genelist, 
                         TERM2GENE = TERM2GENE, 
                         minGSSize = 15, 
                         maxGSSize = 500, 
                         eps = 1e-100,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = 'BH',
                         by = 'fgsea',
                         seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_Ncoa1_FDR5.pdf"))
gseaplot2(gseaResult.Ncoa1, geneSetID = c("Ncoa1"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()

# Usp7

TERM2GENE <- data.frame(TERM = rep("Usp7", length(my.Usp7.genes)), 
                        GENE = my.Usp7.genes, 
                        stringsAsFactors=FALSE)

gseaResult.Usp7 <- GSEA(MM.FMT.ovary.genelist, 
                        TERM2GENE = TERM2GENE, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        eps = 1e-100,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        by = 'fgsea',
                        seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_Usp7_FDR5.pdf"))
gseaplot2(gseaResult.Usp7, geneSetID = c("Usp7"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()

################################################################################
sink(file = paste(Sys.Date(), "MM_FMT_TFs_GSEA_session_Info.txt", sep ="_"))
sessionInfo()
sink()
