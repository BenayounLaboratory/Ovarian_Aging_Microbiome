options(stringsAsFactors = FALSE)

library('dplyr')

###############################################
# Menopause-microbiome project
# Public ovarian aging bulk mRNA-seq analysis - CRA003645 and PRJNA1002222
# Filter ovarian aging signifiant genes for GSEA
###############################################

###############################################
# 1. Import results from DEseq2
###############################################

# CRA003645 - 12M/3M
my.CRA003645 <- read.table("../Ovarian_aging_CRA003645_preprocessing/DESeq2/2024-04-02_MM_CRA003645_DESeq2_aging_Age_all_genes_statistics.txt")

dim(my.CRA003645)     # [1] 10499     6

# PRJNA1002222 - 12M/2M
my.PRJNA1002222 <- read.table("../Ovarian_aging_PRJNA1002222_preprocessing/DESeq2/2024-04-05_MM_PRJNA1002222_DESeq2_aging_Age_all_genes_statistics.txt")

dim(my.PRJNA1002222)     # [1] 20430     6

###############################################
# 2. Combine results based using Fisher's method
###############################################

# Find common genes between two datasets
common.genes <- intersect(rownames(my.CRA003645), rownames(my.PRJNA1002222))

length(common.genes)     # [1] 8764

# Filter and check order of rownames based on common gene list

my.CRA003645.cl <- my.CRA003645[common.genes,]
my.PRJNA1002222.cl <- my.PRJNA1002222[common.genes,]

identical(rownames(my.CRA003645.cl), rownames(my.PRJNA1002222.cl))     # TRUE

my.CRA003645.cl$genes <- rownames(my.CRA003645.cl)
my.PRJNA1002222.cl$genes <- rownames(my.PRJNA1002222.cl)

merged.dataset <- merge(my.CRA003645.cl, my.PRJNA1002222.cl, by = "genes", suffixes = c(".CRA", ".PRJ"))

# Split dataset based on positive and  negative log2FC

# Positive - Up with ovarian aging

upwithage <- merged.dataset[merged.dataset$log2FoldChange.CRA > 0 & merged.dataset$log2FoldChange.PRJ > 0,]     # 2413 genes

# Negative - Down with ovarian aging

downwithage <- merged.dataset[merged.dataset$log2FoldChange.CRA < 0 & merged.dataset$log2FoldChange.PRJ < 0,]   # 2221 genes


# Calculate p-value + BH correction

upwithage$pvalue <- upwithage$pvalue.CRA * upwithage$pvalue.PRJ
upwithage$pvalue_adj_BH <- p.adjust(upwithage$pvalue, method = "BH")

downwithage$pvalue <- downwithage$pvalue.CRA * downwithage$pvalue.PRJ
downwithage$pvalue_adj_BH <- p.adjust(downwithage$pvalue, method = "BH")

# Export table

my.outprefix <- paste(Sys.Date(),"MM_CRA003645_and_PRJNA1002222_combined_statistics_aging",sep="_")

my.out.stats.all <- paste(my.outprefix, "All_genes_statistics.txt", sep = "_")
write.table(merged.dataset, file = my.out.stats.all, sep = "\t", row.names = F, quote = F)

my.out.stats.Age.up <- paste(my.outprefix,"Up_with_age_all_genes_statistics.txt",sep = "_")
write.table(upwithage, file = my.out.stats.Age.up , sep = "\t" , row.names = F, quote=F)

my.out.stats.Age.down <- paste(my.outprefix,"Down_with_age_all_genes_statistics.txt",sep = "_")
write.table(downwithage, file = my.out.stats.Age.down , sep = "\t" , row.names = F, quote=F)

################################################################################
sink(file = paste(Sys.Date(),"MM_Ovarian_aging_public_data_gene_filtering_session_Info.txt", sep =""))
sessionInfo()
sink()

