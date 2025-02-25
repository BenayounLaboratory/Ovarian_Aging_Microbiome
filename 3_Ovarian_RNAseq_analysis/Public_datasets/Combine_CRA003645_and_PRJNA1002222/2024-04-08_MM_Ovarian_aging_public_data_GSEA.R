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
# Public ovarian aging bulk mRNA-seq analysis - CRA003645 & PRJNA1002222
# Use DGEs for GSEA - FMT cohort
###############################################

# Load ovarian aging DESeq2 data

public.data.statistics.upwithage <- read.table("./2024-04-09_MM_CRA003645_and_PRJNA1002222_combined_statistics_aging_Up_with_age_all_genes_statistics.txt", header = TRUE)
public.data.statistics.downwithage <- read.table("./2024-04-09_MM_CRA003645_and_PRJNA1002222_combined_statistics_aging_Down_with_age_all_genes_statistics.txt", header = TRUE)

# Filter sig genes

upwithage.sig <- public.data.statistics.upwithage$genes[public.data.statistics.upwithage$pvalue_adj_BH < 0.00001]
downwithage.sig <- public.data.statistics.downwithage$genes[public.data.statistics.downwithage$pvalue_adj_BH < 0.00001]

# Load FMT DESeq2 data
load("../../FMT_cohort/DESeq2/2024-03-29_MM_FMT_Ovaries_FMT.RData")

##### Prepare gene list using DESeq2 t-statistic to rank genes #####

my.MM.FMT.ovary <- data.frame(res.FMT)
my.MM.FMT.ovary$Gene.name <- rownames(my.MM.FMT.ovary)

MM.FMT.ovary.genelist <- my.MM.FMT.ovary$stat
names(MM.FMT.ovary.genelist) <- my.MM.FMT.ovary$Gene.name
MM.FMT.ovary.genelist <- sort(MM.FMT.ovary.genelist, decreasing = TRUE)

# Age.UP

TERM2GENE <- data.frame(TERM = rep("Age.UP", length(upwithage.sig)), 
                        GENE = upwithage.sig, 
                        stringsAsFactors=FALSE)

gseaResult.age.UP <- GSEA(MM.FMT.ovary.genelist, 
                          TERM2GENE = TERM2GENE, 
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          eps = 1e-100,
                          pvalueCutoff = 1,
                          pAdjustMethod = 'BH',
                          by = 'fgsea',
                          seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_age_UP.pdf"))
gseaplot2(gseaResult.age.UP, geneSetID = c("Age.UP"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()

# Age.DOWN

TERM2GENE <- data.frame(TERM = rep("Age.DOWN", length(downwithage.sig)), 
                        GENE = downwithage.sig, 
                        stringsAsFactors=FALSE)

gseaResult.age.DOWN <- GSEA(MM.FMT.ovary.genelist, 
                            TERM2GENE = TERM2GENE, 
                            minGSSize = 15, 
                            maxGSSize = 500, 
                            eps = 1e-100,
                            pvalueCutoff = 1,
                            pAdjustMethod = 'BH',
                            by = 'fgsea',
                            seed = TRUE)

pdf(paste(Sys.Date(),"_MM_FMT_RNAseq_GSEA_age_DOWN.pdf"))
gseaplot2(gseaResult.age.DOWN, geneSetID = c("Age.DOWN"), color = c('red4'), base_size = 9, pvalue_table = TRUE)
dev.off()

################################################################################
sink(file = paste(Sys.Date(), "MM_FMT_GSEA_session_Info.txt", sep ="_"))
sessionInfo()
sink()
