setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Revision_nat_aging/Data/VCD_combined/4_Alpha_and_beta_diversity_analysis/")

library('qiime2R')
library(tidyverse)
library(mixOmics)
library(PLSDAbatch)

##############################################
# Menopause-microbiome project - revision
# VCD cohort
# Use CLR transformation for pre-processing
# https://evayiwenwang.github.io/Managing_batch_effects/index.html#data-processing
##############################################

##############################################
# 1. Import data files
##############################################

# Import files

metadata.VCD <- read_q2metadata("../1_Manifest_and_metadata/MM_VCD_cohort_combined_sample-metadata.tsv")
metadata.VCD$treatment <- factor(metadata.VCD$treatment, levels = c("CTL", "VCD"))

SVs.VCD      <- read_qza("~/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Revision_nat_aging/Data/VCD_combined/3_Feature_table_and_rep_seq_table/MM_VCD_combined_table.qza")$data

###############################
# 2. Pre-processing
###############################

# 1. Prefilter the count data to remove features with excess zeroes across all samples

VCD.index.keep <- which(colSums(SVs.VCD)*100/(sum(colSums(SVs.VCD))) > 0.01)
VCD.count.keep <- t(SVs.VCD[, VCD.index.keep])
dim(VCD.count.keep)
# [1]   21 1684

# 2. Add offset of 1 for CLR transformation

SVs.VCD <- VCD.count.keep +1

# 3. Centered log-ratio transformation

MM.VCD.clr <- logratio.transfo(SVs.VCD, logratio = 'CLR')
class(MM.VCD.clr) <- 'matrix' 

###############################
# 3. Plot PCA
###############################

VCD.pca <- pca(MM.VCD.clr, ncomp = 3)

pdf(paste(Sys.Date(),"MM_VCD_combined_PCA_CLR.pdf", sep = "_"), height = 10, width = 10)
Scatter_Density(object = VCD.pca, batch = metadata.VCD$sex, trt = metadata.VCD$treatment,
                trt.legend.title = 'VCD', 
                title = 'MM VCD PCA - after CLR')
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_VCD_CLR_PCA_plot_generation.txt", sep =""))
sessionInfo()
sink()
