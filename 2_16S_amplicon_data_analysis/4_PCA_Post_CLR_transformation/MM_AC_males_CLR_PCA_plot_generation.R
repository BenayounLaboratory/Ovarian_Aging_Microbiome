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
# AC cohort - males
# Use CLR transformation for pre-processing
# https://evayiwenwang.github.io/Managing_batch_effects/index.html#data-processing
##############################################

##############################################
# 1. Import data files
##############################################

# Import files

MM.males.AC16.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC16_males_table.qza")
MM.males.AC17.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC17_males_table.qza")
MM.males.AC22.table <- read_qza("../3_Feature_table_and_rep_seq_table/MM_AC22_males_table.qza")

MM.males.AC16.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC16_AC_males_sample-metadata.txt")
MM.males.AC17.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC17_AC_males_sample-metadata.txt")
MM.males.AC22.metadata <- read_q2metadata("../1_Manifest_and_metadata/MM_AC22_AC_males_sample-metadata.txt")

MM.males.AC16.count <- MM.males.AC16.table$data
MM.males.AC17.count <- MM.males.AC17.table$data
MM.males.AC22.count <- MM.males.AC22.table$data

# Combine data from cohorts AC 16, 17 and 22

# Find common features
common.features <- Reduce(intersect, list(rownames(MM.males.AC16.count), rownames(MM.males.AC17.count), rownames(MM.males.AC22.count)))

# Filter tables to only include common features
AC16.common <- MM.males.AC16.count[common.features,]
AC17.common <- MM.males.AC17.count[common.features,]
AC22.common <- MM.males.AC22.count[common.features,]

# Check that order of the row names are consistent

identical(rownames(AC16.common), rownames(AC17.common))     # TRUE
identical(rownames(AC16.common), rownames(AC22.common))     # TRUE

# Combine tables
MM.males.AC.combined.counts <- cbind(AC16.common, AC17.common, AC22.common)

MM.males.AC.combined.counts <- t(MM.males.AC.combined.counts)
MM.males.AC.combined.counts <- as.data.frame(MM.males.AC.combined.counts)

# Combine metadata

MM.males.AC.combined.metadata <- rbind(MM.males.AC16.metadata, MM.males.AC17.metadata, MM.males.AC22.metadata)
rownames(MM.males.AC.combined.metadata) <- MM.males.AC.combined.metadata$SampleID

# Re-order count and metadata to match rowname orders

matched.rownames <- rownames(MM.males.AC.combined.metadata)
MM.males.AC.combined.counts <- MM.males.AC.combined.counts[match(matched.rownames, rownames(MM.males.AC.combined.counts)), , drop = FALSE]

identical(rownames(MM.males.AC.combined.counts), rownames(MM.males.AC.combined.metadata))     # TRUE

###############################
# 2. Pre-processing
###############################

# 1. Prefilter the count data to remove features with excess zeroes across all samples

AC.index.keep <- which(colSums(MM.males.AC.combined.counts)*100/(sum(colSums(MM.males.AC.combined.counts))) > 0.01)
AC.count.keep <- MM.males.AC.combined.counts[, AC.index.keep]
dim(AC.count.keep)
# [1]  28 315

# 2. Add offset of 1 for CLR transformation

AC.count.keep <- AC.count.keep +1

# 3. Centered log-ratio transformation

MM.males.AC.clr <- logratio.transfo(AC.count.keep, logratio = 'CLR')
class(MM.males.AC.clr) <- 'matrix' 

###############################
# 3. Batch detection
###############################

AC.pca.before <- pca(MM.males.AC.clr, ncomp = 3)

Scatter_Density(object = AC.pca.before, batch = MM.males.AC.combined.metadata$batch, 
                trt = MM.males.AC.combined.metadata$age_cat, 
                batch.legend.title = 'Batch', 
                trt.legend.title = 'Age', 
                title = 'Before batch effect correction (AC#16,17,22)')

###############################
# 4. Test batch correction methods
###############################

# limma

MM.males.AC.combined.mod <- model.matrix( ~ age_cat, MM.males.AC.combined.metadata) # full model
MM.males.AC.combined.mod0 <- model.matrix( ~ 1, data = MM.males.AC.combined.metadata) # null model
MM.males.AC.combined.sva.n <- num.sv(dat = t(MM.males.AC.clr), mod = MM.males.AC.combined.mod, method = "b")     #2

MM.males.AC.combined.sva <- sva(dat = t(MM.males.AC.clr), mod = MM.males.AC.combined.mod, 
                                mod0 = MM.males.AC.combined.mod0, n.sv = MM.males.AC.combined.sva.n)

MM.males.AC.combined.mod.bat <- cbind(MM.males.AC.combined.mod, MM.males.AC.combined.sva$sv)
MM.males.AC.combined.mod0.bat <- cbind(MM.males.AC.combined.mod0, MM.males.AC.combined.sva$sv)

MM.males.AC.combined.sva.trt_p <- f.pvalue(t(MM.males.AC.clr), MM.males.AC.combined.mod.bat, MM.males.AC.combined.mod0.bat)
MM.males.AC.combined.sva.trt_adjp <- p.adjust(MM.males.AC.combined.sva.trt_p, method='fdr')

MM.males.AC.combined.limma <- t(removeBatchEffect(t(MM.males.AC.clr), batch = MM.males.AC.combined.metadata$batch, 
                                                  design = MM.males.AC.combined.mod))

# Combat

MM.males.AC.combined.combat <- t(ComBat(t(MM.males.AC.clr), batch = MM.males.AC.combined.metadata$batch, 
                                        mod = MM.males.AC.combined.mod, par.prior = F, prior.plots = F))

# FAbatch

# sponge data
MM.males.AC.combined.fabatch.obj <- fabatch(x = MM.males.AC.clr, 
                                            y = as.factor(as.numeric(MM.males.AC.combined.metadata$age_num)), 
                                            batch = MM.males.AC.combined.metadata$batch)

MM.males.AC.combined.fabatch <- MM.males.AC.combined.fabatch.obj$xadj


###############################
## plot results
###############################

MM.males.AC.combined.pca.before <- pca(MM.males.AC.clr, ncomp = 3)
MM.males.AC.combined.pca.limma <- pca(MM.males.AC.combined.limma, ncomp = 3)
MM.males.AC.combined.pca.combat <- pca(MM.males.AC.combined.combat, ncomp = 3)
MM.males.AC.combined.pca.fabatch <- pca(MM.males.AC.combined.fabatch, ncomp = 3)


p.MM.males.AC.combined.pca.before <- Scatter_Density(MM.males.AC.combined.pca.before, batch = MM.males.AC.combined.metadata$batch, trt = MM.males.AC.combined.metadata$age_cat,
                                                     batch.legend.title = 'Cohort (batch)', 
                                                     trt.legend.title = 'Age (trt)', 
                                                     title = 'Before batch effect correction (Cohort_combined_AC)',
                                                     title.cex = 0.7,
                                                     legend.cex = 0.5,
                                                     legend.title.cex = 0.5)

p.MM.males.AC.combined.pca.limma <- Scatter_Density(MM.males.AC.combined.pca.limma, batch = MM.males.AC.combined.metadata$batch, trt = MM.males.AC.combined.metadata$age_cat,
                                                    batch.legend.title = 'Cohort (batch)', 
                                                    trt.legend.title = 'Age (trt)', 
                                                    title = 'Batch correction - limma (Cohort_combined_AC)',
                                                    title.cex = 0.7,
                                                    legend.cex = 0.5,
                                                    legend.title.cex = 0.5)

p.MM.males.AC.combined.pca.combat <- Scatter_Density(MM.males.AC.combined.pca.combat, batch = MM.males.AC.combined.metadata$batch, trt = MM.males.AC.combined.metadata$age_cat,
                                                     batch.legend.title = 'Cohort (batch)', 
                                                     trt.legend.title = 'Age (trt)', 
                                                     title = 'Batch correction - combat (Cohort_combined_AC)',
                                                     title.cex = 0.7,
                                                     legend.cex = 0.5,
                                                     legend.title.cex = 0.5)

p.MM.males.AC.combined.pca.fabatch <- Scatter_Density(MM.males.AC.combined.pca.fabatch, batch = MM.males.AC.combined.metadata$batch, trt = MM.males.AC.combined.metadata$age_cat,
                                                      batch.legend.title = 'Cohort (batch)', 
                                                      trt.legend.title = 'Age (trt)', 
                                                      title = 'Batch correction - FAbatch (Cohort_combined_AC)',
                                                      title.cex = 0.7,
                                                      legend.cex = 0.5,
                                                      legend.title.cex = 0.5)

pdf(paste(Sys.Date(),"MM_AC_males_combined_AC_PCA_batch_correction.pdf", sep = "_"), height = 7, width = 10)
grid.arrange(p.MM.males.AC.combined.pca.before, p.MM.males.AC.combined.pca.limma, p.MM.males.AC.combined.pca.combat, p.MM.males.AC.combined.pca.fabatch, ncol = 2)
dev.off()

###############################
## Variance calculation
###############################

# MM.males.AC data
MM.males.AC.form <- ~ MM.males.AC.combined.metadata$age_cat + MM.males.AC.combined.metadata$batch
MM.males.AC.info <- as.data.frame(cbind(rownames(MM.males.AC.clr), MM.males.AC.combined.metadata$age_cat,  MM.males.AC.combined.metadata$batch))
rownames(MM.males.AC.info) <- rownames(MM.males.AC.clr)

# before
MM.males.AC.varPart.before <- fitExtractVarPartModel(exprObj = t(MM.males.AC.clr), 
                                                     formula = MM.males.AC.form, 
                                                     data = MM.males.AC.info)

# removeBatchEffect
MM.males.AC.varPart.limma <- fitExtractVarPartModel(exprObj = t(MM.males.AC.combined.limma), 
                                                    formula = MM.males.AC.form, 
                                                    data = MM.males.AC.info)

# combat
MM.males.AC.varPart.combat <- fitExtractVarPartModel(exprObj = t(MM.males.AC.combined.combat), 
                                                     formula = MM.males.AC.form, 
                                                     data = MM.males.AC.info)

# FAbatch
MM.males.AC.varPart.fabatch <- fitExtractVarPartModel(exprObj = t(MM.males.AC.combined.fabatch), 
                                                      formula = MM.males.AC.form, 
                                                      data = MM.males.AC.info)

# Extract the variance of trt and batch

# before
MM.males.AC.varmat.before <- as.matrix(MM.males.AC.varPart.before[ ,1:2])

# removeBatchEffect
MM.males.AC.varmat.limma <- as.matrix(MM.males.AC.varPart.limma[ ,1:2])

# ComBat
MM.males.AC.varmat.combat <- as.matrix(MM.males.AC.varPart.combat[ ,1:2])

# FAbatch
MM.males.AC.varmat.fabatch <- as.matrix(MM.males.AC.varPart.fabatch[ ,1:2])

# merge results
MM.males.AC.variance <- c(as.vector(MM.males.AC.varmat.before), as.vector(MM.males.AC.varmat.limma),
                          as.vector(MM.males.AC.varmat.combat), as.vector(MM.males.AC.varmat.fabatch))

# add batch, trt and methods info
MM.males.AC.variance <- cbind(variance = MM.males.AC.variance, 
                              Type = rep(c('Tissue', 'Batch'), each = ncol(MM.males.AC.clr)),
                              method = rep(c('Before', 'limma', 'ComBat', 'FAbatch'), 
                                           each = 2*ncol(MM.males.AC.clr)))

# reorder levels  
MM.males.AC.variance <- as.data.frame(MM.males.AC.variance)
MM.males.AC.variance$method <- factor(MM.males.AC.variance$method, 
                                      levels = unique(MM.males.AC.variance$method))
MM.males.AC.variance$variance <- as.numeric(as.character(MM.males.AC.variance$variance))

pdf(paste(Sys.Date(),"MM_AC_cohort_males_combined_AC_PCA_batch_correction_variance_calculation.pdf", sep = "_"), height = 7, width = 10)
ggplot(MM.males.AC.variance, aes(x = Type, y = variance, fill = Type)) + 
  geom_boxplot() + facet_grid(cols = vars(method)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        strip.text = element_text(size = 12), panel.grid = element_blank(), 
        axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) + 
  labs(x = 'Type', y = 'Proportion Variance', name = 'Type') + ylim(0,1)
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_AC_males_CLR_batch_correction_benchmarking.txt", sep =""))
sessionInfo()
sink()
