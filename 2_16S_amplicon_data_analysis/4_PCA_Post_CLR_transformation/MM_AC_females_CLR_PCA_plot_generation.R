library('qiime2R')
library(tidyverse)
library(mixOmics)
library(PLSDAbatch)
library(limma)
library(sva) 
library(bapred) 
library(gridExtra)
library(variancePartition)

##############################################
# Menopause-microbiome project
# AC cohort - females
# Use CLR transformation for pre-processing
# https://evayiwenwang.github.io/Managing_batch_effects/index.html#data-processing
##############################################

##############################################
# 1. Import data files
##############################################

# Import files

MM.females.AC16.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC16_females_table.qza")
MM.females.AC17.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC17_females_table.qza")
MM.females.AC22.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC22_females_table.qza")
MM.females.AC24.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC24_females_table.qza")

MM.females.AC16.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC16_AC_females_sample-metadata.txt")
MM.females.AC17.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC17_AC_females_sample-metadata.txt")
MM.females.AC22.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC22_AC_females_sample-metadata.txt")
MM.females.AC24.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC24_AC_females_sample-metadata.txt")

MM.females.AC16.count <- MM.females.AC16.table$data
MM.females.AC17.count <- MM.females.AC17.table$data
MM.females.AC22.count <- MM.females.AC22.table$data
MM.females.AC24.count <- MM.females.AC24.table$data

# Combine data from cohorts AC 16, 17, 22 and 24

# Find common features
common.features <- Reduce(intersect, list(rownames(MM.females.AC16.count), rownames(MM.females.AC17.count), rownames(MM.females.AC22.count), rownames(MM.females.AC24.count)))

# Filter tables to only include common features
AC16.common <- MM.females.AC16.count[common.features,]
AC17.common <- MM.females.AC17.count[common.features,]
AC22.common <- MM.females.AC22.count[common.features,]
AC24.common <- MM.females.AC24.count[common.features,]

# Check that order of the row names are consistent

identical(rownames(AC16.common), rownames(AC17.common))     # TRUE
identical(rownames(AC16.common), rownames(AC22.common))     # TRUE
identical(rownames(AC16.common), rownames(AC24.common))     # TRUE

# Combine tables
MM.females.AC.combined.counts <- cbind(AC16.common, AC17.common, AC22.common, AC24.common)

MM.females.AC.combined.counts <- t(MM.females.AC.combined.counts)
MM.females.AC.combined.counts <- as.data.frame(MM.females.AC.combined.counts)

# Combine metadata

MM.females.AC.combined.metadata <- rbind(MM.females.AC16.metadata, MM.females.AC17.metadata, MM.females.AC22.metadata, MM.females.AC24.metadata)
rownames(MM.females.AC.combined.metadata) <- MM.females.AC.combined.metadata$SampleID

# Re-order count and metadata to match rowname orders

matched.rownames <- rownames(MM.females.AC.combined.metadata)
MM.females.AC.combined.counts <- MM.females.AC.combined.counts[match(matched.rownames, rownames(MM.females.AC.combined.counts)), , drop = FALSE]

identical(rownames(MM.females.AC.combined.counts), rownames(MM.females.AC.combined.metadata))     # TRUE

###############################
# 2. Pre-processing
###############################

# 1. Prefilter the count data to remove features with excess zeroes across all samples

AC.index.keep <- which(colSums(MM.females.AC.combined.counts)*100/(sum(colSums(MM.females.AC.combined.counts))) > 0.01)
AC.count.keep <- MM.females.AC.combined.counts[, AC.index.keep]
dim(AC.count.keep)
# [1]  39 266

# 2. Add offset of 1 for CLR transformation

AC.count.keep <- AC.count.keep +1

# 3. Centered log-ratio transformation

MM.females.AC.clr <- logratio.transfo(AC.count.keep, logratio = 'CLR')
class(MM.females.AC.clr) <- 'matrix' 

###############################
# 3. Batch detection
###############################

AC.pca.before <- pca(MM.females.AC.clr, ncomp = 3)

Scatter_Density(object = AC.pca.before, batch = MM.females.AC.combined.metadata$batch, 
                trt = MM.females.AC.combined.metadata$age_cat, 
                batch.legend.title = 'Batch', 
                trt.legend.title = 'Age', 
                title = 'Before batch effect correction (AC#16,17,22,24)')

###############################
# 4. Test batch correction methods
###############################

# limma

MM.females.AC.combined.mod <- model.matrix( ~ age_cat, MM.females.AC.combined.metadata) # full model
MM.females.AC.combined.mod0 <- model.matrix( ~ 1, data = MM.females.AC.combined.metadata) # null model
MM.females.AC.combined.sva.n <- num.sv(dat = t(MM.females.AC.clr), mod = MM.females.AC.combined.mod, method = "b")     #4

MM.females.AC.combined.sva <- sva(dat = t(MM.females.AC.clr), mod = MM.females.AC.combined.mod, 
                                  mod0 = MM.females.AC.combined.mod0, n.sv = MM.females.AC.combined.sva.n)

MM.females.AC.combined.mod.bat <- cbind(MM.females.AC.combined.mod, MM.females.AC.combined.sva$sv)
MM.females.AC.combined.mod0.bat <- cbind(MM.females.AC.combined.mod0, MM.females.AC.combined.sva$sv)

MM.females.AC.combined.sva.trt_p <- f.pvalue(t(MM.females.AC.clr), MM.females.AC.combined.mod.bat, MM.females.AC.combined.mod0.bat)
MM.females.AC.combined.sva.trt_adjp <- p.adjust(MM.females.AC.combined.sva.trt_p, method='fdr')

MM.females.AC.combined.limma <- t(removeBatchEffect(t(MM.females.AC.clr), batch = MM.females.AC.combined.metadata$batch, 
                                                    design = MM.females.AC.combined.mod))

# Combat

MM.females.AC.combined.combat <- t(ComBat(t(MM.females.AC.clr), batch = MM.females.AC.combined.metadata$batch, 
                                          mod = MM.females.AC.combined.mod, par.prior = F, prior.plots = F))

# FAbatch

# sponge data
MM.females.AC.combined.fabatch.obj <- fabatch(x = MM.females.AC.clr, 
                                              y = as.factor(as.numeric(MM.females.AC.combined.metadata$age_num)), 
                                              batch = MM.females.AC.combined.metadata$batch)

MM.females.AC.combined.fabatch <- MM.females.AC.combined.fabatch.obj$xadj


###############################
## plot results
###############################

MM.females.AC.combined.pca.before <- pca(MM.females.AC.clr, ncomp = 3)
MM.females.AC.combined.pca.limma <- pca(MM.females.AC.combined.limma, ncomp = 3)
MM.females.AC.combined.pca.combat <- pca(MM.females.AC.combined.combat, ncomp = 3)
MM.females.AC.combined.pca.fabatch <- pca(MM.females.AC.combined.fabatch, ncomp = 3)

p.MM.females.AC.combined.pca.before <- Scatter_Density(MM.females.AC.combined.pca.before, batch = MM.females.AC.combined.metadata$batch, trt = MM.females.AC.combined.metadata$age_cat,
                                                       batch.legend.title = 'Cohort (batch)', 
                                                       trt.legend.title = 'Age (trt)', 
                                                       title = 'Before batch effect correction (Cohort_combined_AC)',
                                                       title.cex = 0.7,
                                                       legend.cex = 0.5,
                                                       legend.title.cex = 0.5)

p.MM.females.AC.combined.pca.limma <- Scatter_Density(MM.females.AC.combined.pca.limma, batch = MM.females.AC.combined.metadata$batch, trt = MM.females.AC.combined.metadata$age_cat,
                                                      batch.legend.title = 'Cohort (batch)', 
                                                      trt.legend.title = 'Age (trt)', 
                                                      title = 'Batch correction - limma (Cohort_combined_AC)',
                                                      title.cex = 0.7,
                                                      legend.cex = 0.5,
                                                      legend.title.cex = 0.5)

p.MM.females.AC.combined.pca.combat <- Scatter_Density(MM.females.AC.combined.pca.combat, batch = MM.females.AC.combined.metadata$batch, trt = MM.females.AC.combined.metadata$age_cat,
                                                       batch.legend.title = 'Cohort (batch)', 
                                                       trt.legend.title = 'Age (trt)', 
                                                       title = 'Batch correction - combat (Cohort_combined_AC)',
                                                       title.cex = 0.7,
                                                       legend.cex = 0.5,
                                                       legend.title.cex = 0.5)

p.MM.females.AC.combined.pca.fabatch <- Scatter_Density(MM.females.AC.combined.pca.fabatch, batch = MM.females.AC.combined.metadata$batch, trt = MM.females.AC.combined.metadata$age_cat,
                                                        batch.legend.title = 'Cohort (batch)', 
                                                        trt.legend.title = 'Age (trt)', 
                                                        title = 'Batch correction - FAbatch (Cohort_combined_AC)',
                                                        title.cex = 0.7,
                                                        legend.cex = 0.5,
                                                        legend.title.cex = 0.5)

pdf(paste(Sys.Date(),"MM_AC_females_combined_AC_PCA_batch_correction.pdf", sep = "_"), height = 7, width = 10)
grid.arrange(p.MM.females.AC.combined.pca.before, p.MM.females.AC.combined.pca.limma, p.MM.females.AC.combined.pca.combat, p.MM.females.AC.combined.pca.fabatch, ncol = 2)
dev.off()

###############################
## Variance calculation
###############################

# MM.females.AC data
MM.females.AC.form <- ~ MM.females.AC.combined.metadata$age_cat + MM.females.AC.combined.metadata$batch
MM.females.AC.info <- as.data.frame(cbind(rownames(MM.females.AC.clr), MM.females.AC.combined.metadata$age_cat,  MM.females.AC.combined.metadata$batch))
rownames(MM.females.AC.info) <- rownames(MM.females.AC.clr)

# before
MM.females.AC.varPart.before <- fitExtractVarPartModel(exprObj = t(MM.females.AC.clr), 
                                                       formula = MM.females.AC.form, 
                                                       data = MM.females.AC.info)

# removeBatchEffect
MM.females.AC.varPart.limma <- fitExtractVarPartModel(exprObj = t(MM.females.AC.combined.limma), 
                                                      formula = MM.females.AC.form, 
                                                      data = MM.females.AC.info)

# combat
MM.females.AC.varPart.combat <- fitExtractVarPartModel(exprObj = t(MM.females.AC.combined.combat), 
                                                       formula = MM.females.AC.form, 
                                                       data = MM.females.AC.info)

# FAbatch
MM.females.AC.varPart.fabatch <- fitExtractVarPartModel(exprObj = t(MM.females.AC.combined.fabatch), 
                                                        formula = MM.females.AC.form, 
                                                        data = MM.females.AC.info)

# Extract the variance of trt and batch

# before
MM.females.AC.varmat.before <- as.matrix(MM.females.AC.varPart.before[ ,1:2])

# removeBatchEffect
MM.females.AC.varmat.limma <- as.matrix(MM.females.AC.varPart.limma[ ,1:2])

# ComBat
MM.females.AC.varmat.combat <- as.matrix(MM.females.AC.varPart.combat[ ,1:2])

# FAbatch
MM.females.AC.varmat.fabatch <- as.matrix(MM.females.AC.varPart.fabatch[ ,1:2])

# merge results
MM.females.AC.variance <- c(as.vector(MM.females.AC.varmat.before), as.vector(MM.females.AC.varmat.limma),
                            as.vector(MM.females.AC.varmat.combat), as.vector(MM.females.AC.varmat.fabatch))

# add batch, trt and methods info
MM.females.AC.variance <- cbind(variance = MM.females.AC.variance, 
                                Type = rep(c('Tissue', 'Batch'), each = ncol(MM.females.AC.clr)),
                                method = rep(c('Before', 'limma', 'ComBat', 'FAbatch'), 
                                             each = 2*ncol(MM.females.AC.clr)))

# reorder levels  
MM.females.AC.variance <- as.data.frame(MM.females.AC.variance)
MM.females.AC.variance$method <- factor(MM.females.AC.variance$method, 
                                        levels = unique(MM.females.AC.variance$method))
MM.females.AC.variance$variance <- as.numeric(as.character(MM.females.AC.variance$variance))

pdf(paste(Sys.Date(),"MM_AC_cohort_females_combined_AC_PCA_batch_correction_variance_calculation.pdf", sep = "_"), height = 7, width = 10)
ggplot(MM.females.AC.variance, aes(x = Type, y = variance, fill = Type)) + 
  geom_boxplot() + facet_grid(cols = vars(method)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        strip.text = element_text(size = 12), panel.grid = element_blank(), 
        axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) + 
  labs(x = 'Type', y = 'Proportion Variance', name = 'Type') + ylim(0,1)
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_AC_females_CLR_batch_correction_benchmarking.txt", sep =""))
sessionInfo()
sink()
