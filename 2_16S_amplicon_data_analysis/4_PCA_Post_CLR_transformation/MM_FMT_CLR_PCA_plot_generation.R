library('qiime2R')
library(tidyverse)
library(mixOmics)
library(PLSDAbatch)

##############################################
# Menopause-microbiome project
# FMT cohort
# Use CLR transformation for pre-processing
# https://evayiwenwang.github.io/Managing_batch_effects/index.html#data-processing
##############################################

##############################################
# 1. Import data files
##############################################

# Import files

metadata.FMT <- read_q2metadata("../1_Manifest_and_metadata/MM_FMT_AC_PostFMT_sample-metadata.tsv")
metadata.FMT$fmt <- factor(metadata.FMT$fmt, levels = c("FMT_Y", "FMT_O"))

SVs.FMT      <- read_qza("../3_Feature_table_and_rep_seq_table/MM_FMT_AC_PostFMT_table.qza")$data

###############################
# 2. Pre-processing
###############################

# 1. Prefilter the count data to remove features with excess zeroes across all samples

FMT.index.keep <- which(colSums(SVs.FMT)*100/(sum(colSums(SVs.FMT))) > 0.01)
FMT.count.keep <- t(SVs.FMT[, FMT.index.keep])
dim(FMT.count.keep)
# [1]   14 1511

# 2. Add offset of 1 for CLR transformation

SVs.FMT <- FMT.count.keep +1

# 3. Centered log-ratio transformation

MM.FMT.clr <- logratio.transfo(SVs.FMT, logratio = 'CLR')
class(MM.FMT.clr) <- 'matrix' 

###############################
# 3. Plot PCA
###############################

FMT.pca <- pca(MM.FMT.clr, ncomp = 3)

pdf(paste(Sys.Date(),"MM_FMT_PCA_CLR.pdf", sep = "_"), height = 7, width = 10)
Scatter_Density(object = FMT.pca, trt = metadata.FMT$fmt,
                trt.legend.title = 'FMT', 
                title = 'MM FMT PCA - after CLR')
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_FMT_CLR_PCA_plot.txt", sep =""))
sessionInfo()
sink()
