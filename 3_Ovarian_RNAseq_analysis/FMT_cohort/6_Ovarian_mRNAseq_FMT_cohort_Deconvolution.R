setwd('/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Github/Ovarian_mRNAseq/FMT_cohort/Deconvolution')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library('CSCDRNA')
library('GenomicFeatures')
library('granulator') 
library("Biobase")
library('Seurat')
library('pheatmap')  
library('beeswarm')
library('ggplot2')
library("dplyr")
library(tidyverse)

###############################################
# Menopause-Microbiome project
# Ovarian mRNA-seq - FMT cohort (FMT-YF vs. FMT-EF)
# Run deconvolution - immune vs. non-immune cells
# CSCDRNA & Granulator
###############################################

##### Load MM FMT ovarian RNAseq #####

my.ovary.counts1 <- read.table("../MM_FMT_ovary_mRNAseq_featureCounts_output.txt", sep="\t", header=T)
my.ovary.counts2 <- my.ovary.counts1[,-c(2:6)]

colnames(my.ovary.counts2) <- c("GeneName", 
                                "FMT_YF_1",
                                "FMT_YF_2",
                                "FMT_EF_1",
                                "FMT_EF_2",
                                "FMT_EF_3",
                                "FMT_YF_3",
                                "FMT_YF_4",
                                "FMT_YF_5",
                                "FMT_YF_6",
                                "FMT_YF_7",
                                "FMT_EF_4",
                                "FMT_EF_5",
                                "FMT_EF_6",
                                "FMT_EF_7")

my.ovary.counts <- my.ovary.counts2[,c("GeneName",
                                       "FMT_YF_1",
                                       "FMT_YF_2",
                                       "FMT_YF_3",
                                       "FMT_YF_4",
                                       "FMT_YF_5",
                                       "FMT_YF_6",
                                       "FMT_YF_7",
                                       "FMT_EF_1",
                                       "FMT_EF_2",
                                       "FMT_EF_3",
                                       "FMT_EF_4",
                                       "FMT_EF_5",
                                       "FMT_EF_6",
                                       "FMT_EF_7")]

##### Load reference data - mouse cycling ovary scRNAseq data #####
# https://elifesciences.org/articles/77239
# https://osf.io/924fz/

# Load Seurat object

my.cycling.ovary.scRNAseq <- readRDS("./Murine_cycling_ovarian_scRNAseq_Seurat_object.rds")
my.cycling.ovary.scRNAseq
# An object of class Seurat 
# 23657 features across 34712 samples within 1 assay 
# Active assay: RNA (23657 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

my.cycling.ovary.scRNAseq.counts.table <- as.matrix(GetAssayData(object = my.cycling.ovary.scRNAseq, slot = "counts"))
dim(my.cycling.ovary.scRNAseq.counts.table)
# [1] 23657 34712

###############################################
# Deconvolution using CSCDRNA
###############################################

##### Build ExpressionSet objects with bulk RNAseq and scRNAseq datasets #####

#Build ExpressionSet with bulk data.
MM.FMT.bulk.eset <- ExpressionSet(assayData = as.matrix(my.ovary.counts))

#Build ExpressionSet with single-cell data.
sc.counts.matrix = my.cycling.ovary.scRNAseq.counts.table
individual.labels = my.cycling.ovary.scRNAseq@meta.data$condition
cell.type.labels = my.cycling.ovary.scRNAseq@meta.data$Level0
sample.ids <- colnames(sc.counts.matrix)
#individual.labels and cell.types should be in the same order as in sample.ids.

sc.pheno <- data.frame(check.names = FALSE, check.rows = FALSE,
                       stringsAsFactors = FALSE, row.names = sample.ids,
                       SubjectName = individual.labels, cellType = cell.type.labels)
sc.meta <- data.frame(labelDescription = c("SubjectName", "cellType"),
                      row.names = c("SubjectName", "cellType"))
sc.pdata <- new("AnnotatedDataFrame", data = sc.pheno, varMetadata = sc.meta)
sc.eset <- ExpressionSet(assayData = sc.counts.matrix, phenoData = sc.pdata)

#Run CSCD
analysis <- CSCD(bulk.eset = MM.FMT.bulk.eset, sc.eset = sc.eset,
                 min.p = 0.3, markers = NULL, cell.types = "cellType",
                 subj.names = "SubjectName", verbose=TRUE)

#Estimated cell-type proportions
analysis$bulk.props
#                 FMT_YF_1     FMT_YF_2     FMT_YF_3     FMT_YF_4     FMT_YF_5    FMT_YF_6     FMT_YF_7    FMT_EF_1    FMT_EF_2    FMT_EF_3    FMT_EF_4     FMT_EF_5    FMT_EF_6    FMT_EF_7
# Granulosa   0.5288415759 0.5491150189 0.4993841153 0.5026769196 0.5764349428 0.513143119 0.5233780911 0.477215821 0.477858304 0.503114906 0.478261090 0.4759377768 0.489071543 0.492465264
# Mesenchyme  0.3005745677 0.2938505685 0.3124266080 0.3044897269 0.2888732534 0.304522177 0.3035609719 0.306628356 0.309248106 0.306254638 0.305687839 0.3117687335 0.303103114 0.302293484
# Immune      0.0458338665 0.0430309864 0.0479037731 0.0487783119 0.0402740327 0.047154012 0.0455294854 0.052081066 0.051469454 0.048167857 0.051917821 0.0508130474 0.051026523 0.050952220
# Endothelium 0.0933609670 0.0844859088 0.1049415713 0.1093090132 0.0705922024 0.101476686 0.0964854357 0.125073018 0.122744873 0.107237704 0.126359544 0.1242365413 0.120710457 0.116533603
# Epithelium  0.0310174601 0.0290624272 0.0348721331 0.0338756487 0.0234626933 0.032993661 0.0305284808 0.037818254 0.037676527 0.034463620 0.036509742 0.0363705112 0.034881468 0.036590445
# Oocyte      0.0003715628 0.0004550901 0.0004717994 0.0008703797 0.0003628754 0.000710345 0.0005175351 0.001183485 0.001002737 0.000761275 0.001263964 0.0008733898 0.001206895 0.001164984


# Plot estimated immune cell proportions

immune.cell.prop <- analysis$bulk.props["Immune",]

immune_long <- tibble(Sample = names(immune.cell.prop), 
                      Est.prop = immune.cell.prop) %>%
               mutate(Group = ifelse(grepl("FMT_YF", Sample), "FMT-YF", "FMT-EF"))

immune_long$Group <- factor(immune_long$Group, levels = c("FMT-YF", "FMT-EF"))

# color-code 
my.colors <- rep("#9C509F", 14)
my.colors[grep("FMT-EF",immune_long$Group)] <- "#F99F27"

wilcox.test(immune_long$Est.prop[immune_long$Group == "FMT-YF"], immune_long$Est.prop[immune_long$Group == "FMT-EF"])     # p-value = 0.1166

pdf(paste(Sys.Date(),"_CSCDRNA_deconvolution_pseudobulk_immune_cells_estimated_proportions_FMT-YF_vs_FMT-EF.pdf"), height = 5, width = 6)
beeswarm(Est.prop ~ Group, data = immune_long, pch = 19, pwcol = my.colors,
         ylim = c(0, 0.2))
legend("topright", legend = c("FMT-YF", "FMT-OF"), col = c("#9C509F", "#F99F27"), pch = 19)
dev.off()

###############################################
# Deconvolution using Granulator
###############################################

##### deconvolution analysis accessory functions #####

# Function to generate tpm values
# https://support.bioconductor.org/p/91218/
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# Function to make weighted means of expression
get_wt_mean <- function(exp.matrix,my.weights){
  my.result <- rep(0,dim(exp.matrix)[1])
  for (i in 1: dim(exp.matrix)[1]) {
    my.result[i] <- weighted.mean(exp.matrix[i,], my.weights)
  }
  return(my.result)
}

# Function to generate "fake" tpm values from UMI 
# (length correction of 100 based on 10xGenomics UMI read)
tpmUMI <- function(counts) {
  x <- counts/100
  return(t(t(x)*1e6/colSums(x)))
}

####################################################
# TPM normalization is required by Granulator package
# http://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html

# Generate txdb from gtf using GenomeFeatures
# https://rdrr.io/bioc/GenomicFeatures/man/makeTxDbFromGFF.html
my.GRCm39.txdb <- makeTxDbFromGFF("./GCF_000001635.27_GRCm39_genomic.gff")
my.GRCm39.trscpt.lengths <- transcriptLengths(my.GRCm39.txdb, with.cds_len=T,with.utr5_len=T, with.utr3_len=T)

# Retain only longest transcript for each gene
my.GRCm39.lg.flt <- aggregate(my.GRCm39.trscpt.lengths[,c("nexon","tx_len")], by = list(my.GRCm39.trscpt.lengths$gene_id), FUN = "max")
colnames(my.GRCm39.lg.flt)[1] <- "gene_id"

# Generate a merged matrix with counts and transcript length
my.MM.FMT.ovary.lgth <- merge(my.GRCm39.lg.flt[,c("gene_id","tx_len")], my.ovary.counts, by.x = "gene_id", by.y = "GeneName")
rownames(my.MM.FMT.ovary.lgth) <- my.MM.FMT.ovary.lgth$gene_id

# Generate tpm values
my.MM.FMT.ovary.tpm <- tpm3(my.MM.FMT.ovary.lgth[,-c(1:2)], my.MM.FMT.ovary.lgth$tx_len)

write.table(my.MM.FMT.ovary.tpm, file = paste0(Sys.Date(),"_MM_FMT_ovary_TPM.txt"), sep = "\t", quote = F)

##### Load scRNAseq data (PMID: 36205477) #####
# https://elifesciences.org/articles/77239

######## A. Read in the Seurat and metadata  ########
# Load Seurat object
my.cycling.ovary.scRNAseq
# An object of class Seurat 
# 23657 features across 34712 samples within 1 assay 
# Active assay: RNA (23657 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

table(my.cycling.ovary.scRNAseq@meta.data$Level2)

#           EN_Blood      EN_Blood_Vwf_High  EN_Lymph_Clca3a1_High   EN_Lymph_Clca3a1_Low             Epithelium              GC_Antral             GC_Atretic           GC_CL_Active            GC_CL_Lytic 
#               2938                     84                    197                    282                   1088                   4369                   3080                    670                    635 
#         GC_Estrous             GC_Mitotic               GC_Mural           GC_Preantral        I_Ar_Macrophage               I_B_Cell       I_Dendritic_Cell     I_Gpnmb_Macrophage          I_Granulocyte 
#                879                   2910                    368                   4716                     70                     45                    107                    135                     46 
#     I_Il2ra_T_Cell           I_Macrophage             I_NKT_Cell I_Serpinb10_Macrophage               I_T_Cell     I_Top2a_Macrophage      M_Cortical Stroma  M_Dividing Mesenchyme       M_Immature Theca 
#                 55                    747                     82                     97                    173                     92                   3328                    891                   1880 
# M_Medullary Stroma             M_Pericyte        M_Smooth Muscle  M_Steroidogenic Theca                 Oocyte 
#               1273                    678                   1001                   1774                     22 

# Re-label cell type annotation to generalize cell type
Annotation <- rep("NA", length(my.cycling.ovary.scRNAseq@meta.data$Level2))
Annotation[grep("EN_", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Endothelial"
Annotation[grep("Epithelium", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Epithelial"
Annotation[grep("GC_", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Granulosa"
Annotation[grep("I_", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Immune"
Annotation[grep("M_", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Mesenchymal"
Annotation[grep("Stroma", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Stromal"
Annotation[grep("Theca", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Theca"
Annotation[grep("Oocyte", my.cycling.ovary.scRNAseq@meta.data$Level2)] <- "Oocyte"

Annotation <- data.frame(Annotation)
rownames(Annotation) <- colnames(my.cycling.ovary.scRNAseq@assays$RNA)

# Transfer label to Seurat object
my.cycling.ovary.scRNAseq <- AddMetaData(object = my.cycling.ovary.scRNAseq, metadata = as.vector(Annotation), col.name = "Annotation")

table(my.cycling.ovary.scRNAseq@meta.data$Annotation)
# Endothelial  Epithelial   Granulosa      Immune Mesenchymal      Oocyte     Stromal       Theca 
#        3501        1088       17627        1649        2570          22        4601        3654

# Make global reference pseudobulks
cycling.ovary.pb <- data.frame(AggregateExpression(my.cycling.ovary.scRNAseq, group.by = "Annotation", slot = "counts", assay = "RNA")$RNA)
cycling.ovary.level1.pb <- data.frame(AggregateExpression(my.cycling.ovary.scRNAseq, group.by = "Level1", slot = "counts", assay = "RNA")$RNA)

# Generate tpm values for pure pseudobulk
cycling.ovary.pb.tpm <- tpmUMI(cycling.ovary.pb)
cycling.ovary.level1.pb.tpm <- tpmUMI(cycling.ovary.level1.pb)

##### Make ovary pseudobulks to test deconvolution accuracy #####

# Get "True" Percentages
global.props <- table(my.cycling.ovary.scRNAseq@meta.data$Annotation)/length(my.cycling.ovary.scRNAseq@meta.data$Annotation[!is.na(my.cycling.ovary.scRNAseq@meta.data$Annotation)])
global.level1.props <- table(my.cycling.ovary.scRNAseq@meta.data$Level1)/length(my.cycling.ovary.scRNAseq@meta.data$Level1[!is.na(my.cycling.ovary.scRNAseq@meta.data$Level1)])

#### simulate mixtures
set.seed(123456789)

# Create full pseudobulk for each cell type
cycling.ovary.pure.pb  <- data.frame(AggregateExpression(my.cycling.ovary.scRNAseq, group.by = "Annotation", slot = "counts", assay = "RNA")$RNA)
cycling.ovary.level1.pure.pb  <- data.frame(AggregateExpression(my.cycling.ovary.scRNAseq, group.by = "Level1", slot = "counts", assay = "RNA")$RNA)

# Make up random "shuffles" of ground truth percentages (first 1 is ground truth)
my.pb.freqs <- data.frame("PB_1" = as.numeric(global.props),
                          "PB_2" = sample(as.numeric(global.props)),
                          "PB_3" = sample(as.numeric(global.props)),
                          "PB_4" = sample(as.numeric(global.props)),
                          "PB_5" = sample(as.numeric(global.props)))
rownames(my.pb.freqs) <- names(global.props)

my.level1.pb.freqs <- data.frame("PB_1" = as.numeric(global.level1.props),
                                 "PB_2" = sample(as.numeric(global.level1.props)),
                                 "PB_3" = sample(as.numeric(global.level1.props)),
                                 "PB_4" = sample(as.numeric(global.level1.props)),
                                 "PB_5" = sample(as.numeric(global.level1.props)))
rownames(my.level1.pb.freqs) <- names(global.level1.props)

# Make ovarian bulks using "real" proportions for a weighted mean
my.ov.sim.1 <- get_wt_mean(cycling.ovary.pure.pb, my.pb.freqs$PB_1)
my.ov.sim.2 <- get_wt_mean(cycling.ovary.pure.pb, my.pb.freqs$PB_2)
my.ov.sim.3 <- get_wt_mean(cycling.ovary.pure.pb, my.pb.freqs$PB_3)
my.ov.sim.4 <- get_wt_mean(cycling.ovary.pure.pb, my.pb.freqs$PB_4)
my.ov.sim.5 <- get_wt_mean(cycling.ovary.pure.pb, my.pb.freqs$PB_5)

my.ov.level1.sim.1 <- get_wt_mean(cycling.ovary.level1.pure.pb, my.level1.pb.freqs$PB_1)
my.ov.level1.sim.2 <- get_wt_mean(cycling.ovary.level1.pure.pb, my.level1.pb.freqs$PB_2)
my.ov.level1.sim.3 <- get_wt_mean(cycling.ovary.level1.pure.pb, my.level1.pb.freqs$PB_3)
my.ov.level1.sim.4 <- get_wt_mean(cycling.ovary.level1.pure.pb, my.level1.pb.freqs$PB_4)
my.ov.level1.sim.5 <- get_wt_mean(cycling.ovary.level1.pure.pb, my.level1.pb.freqs$PB_5)

# Make dataframe
cycling.ov.PBmix <- round(data.frame('PB_1'= my.ov.sim.1, 
                                     'PB_2'= my.ov.sim.2, 
                                     'PB_3'= my.ov.sim.3,
                                     'PB_4'= my.ov.sim.4,
                                     'PB_5'= my.ov.sim.5))
rownames(cycling.ov.PBmix) <- rownames(cycling.ovary.pure.pb)

cycling.ov.level1.PBmix <- round(data.frame('PB_1'= my.ov.level1.sim.1, 
                                            'PB_2'= my.ov.level1.sim.2, 
                                            'PB_3'= my.ov.level1.sim.3,
                                            'PB_4'= my.ov.level1.sim.4,
                                            'PB_5'= my.ov.level1.sim.5))
rownames(cycling.ov.level1.PBmix) <- rownames(cycling.ovary.level1.pure.pb)

# Generate tpm values for fake bulks
cycling.ov.PBmix.tpm <- tpmUMI(cycling.ov.PBmix)
cycling.ov.level1.PBmix.tpm <- tpmUMI(cycling.ov.level1.PBmix)

##### Deconvolution #####
# http://bioconductor.org/packages/release/bioc/vignettes/granulator/inst/doc/granulator.html#bulk-pbmcs-rna-seq

# Deconvolution of pseudo bulk RNA-seq data from mouse cycling ovary to determine performance
PBmix.decon <- deconvolute(m = as.matrix(cycling.ov.PBmix), sigMatrix = as.matrix(cycling.ovary.pure.pb))
# WARNING: reaching max number of iterations
PBmix.level1.decon <- deconvolute(m = as.matrix(cycling.ov.level1.PBmix), sigMatrix = as.matrix(cycling.ovary.level1.pure.pb))

save(PBmix.decon, file = "./Granulator_decon_PBmix_cycling_ref_level0.RData")
save(PBmix.level1.decon, file = "./Granulator_decon_PBmix_cycling_ref_level1.RData")

# We can plot the estimated cell type proportions with the function plot_proportions().
# Notice that while the sum of cell types proportions cannot exceed 100%, for some methods part of the bulk RNA-seq signal remains unassigned.
# plot cell type proportions for svr model on ABIS_S0 reference profile
plot_proportions(deconvoluted = PBmix.decon, method = 'svr')
plot_proportions(deconvoluted = PBmix.level1.decon, method = 'svr')

PBmix.bench <- benchmark(deconvoluted = PBmix.decon, ground_truth = as.matrix(t(my.pb.freqs)) )
PBmix.level1.bench <- benchmark(deconvoluted = PBmix.level1.decon, ground_truth = as.matrix(t(my.level1.pb.freqs)) )

# print metrics
PBmix.bench$rank
#   signature  method mean_pcc mean_ccc mean_adj.r2 mean_rmse
# 1      sig1    nnls   0.9991   0.0081      0.9988    0.0018
# 2      sig1     ols   0.9991   0.0081      0.9988    0.0018
# 3      sig1   qprog   0.9991   0.0081      0.9988    0.0018
# 4      sig1 qprogwc   0.9991   0.0081      0.9988    0.0018
# 5      sig1     rls   0.9991   0.0081      0.9988    0.0018
# 6      sig1 dtangle   0.9430   0.0076      0.8638    0.0103
# 7      sig1     svr       NA       NA         NaN    0.0254

PBmix.level1.bench$rank
#   signature  method mean_pcc mean_ccc mean_adj.r2 mean_rmse
# 1      sig1    nnls   0.9966   0.0081      0.9917    0.0020
# 2      sig1     ols   0.9966   0.0081      0.9917    0.0020
# 3      sig1   qprog   0.9966   0.0081      0.9917    0.0020
# 4      sig1 qprogwc   0.9966   0.0081      0.9917    0.0020
# 5      sig1     rls   0.9966   0.0081      0.9917    0.0020
# 6      sig1 dtangle       NA       NA      0.8383    0.0033
# 7      sig1     svr       NA       NA         NaN    0.0133

pdf(paste0(Sys.Date(),"_dothcart_granulator_algorithms_performance_on_PBmix_benchmarking.pdf"), height = 4, width = 4)
dotchart(PBmix.bench$rank$mean_pcc, labels = PBmix.bench$rank$method, xlim = c(0.5,1), pch = 16, xlab = "Mean PCC")
dev.off()

# plot pearson correlation between predictions and true proportions
pdf(paste0(Sys.Date(),"_Deconvolution_algorithms_performance_on_PBmix_PCC_metric.pdf"), height = 5, width = 6)
plot_benchmark(benchmarked = PBmix.bench, metric = 'pcc')
dev.off()

# plot pearson correlation between predictions and true proportions
pdf(paste0(Sys.Date(),"_Deconvolution_algorithms_performance_on_PBmix_level1_annotation_PCC_metric.pdf"), height = 5, width = 6)
plot_benchmark(benchmarked = PBmix.level1.bench, metric = 'pcc')
dev.off()

# Extract regression/correlation plots
pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_NNLS.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'nnls')
dev.off()

pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_OLS.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'ols')
dev.off()

pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_QPROG.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'qprog')
dev.off()

pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_QPROGWC.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'qprogwc')
dev.off()

pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_RLS.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'rls')
dev.off()

pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_DTANGLE.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'dtangle')
dev.off()

pdf(paste0(Sys.Date(),"_realvspred_cell_type_ovary_correlation_SVR.pdf"), height = 5, width = 8)
plot_regress(benchmarked = PBmix.bench, method = 'svr')
dev.off()

##############################################################
##### Run deconvolution on MM FMT ovary data #####

# Deconvolution of pseudo bulk RNA-seq data using all algorithms from benchmark
MM.FMT.decon <- deconvolute(m = as.matrix(my.MM.FMT.ovary.tpm), sigMatrix = as.matrix(cycling.ovary.pure.pb), c("nnls", "ols", "qprog", "qprogwc", "rls", "dtangle", "svr"))
MM.FMT.level1.decon <- deconvolute(m = as.matrix(my.MM.FMT.ovary.tpm), sigMatrix = as.matrix(cycling.ovary.level1.pure.pb), c("nnls", "ols", "qprog", "qprogwc", "rls", "dtangle", "svr"))

save(MM.FMT.decon, file = "./Granulator_decon_MM_FMT_cycling_ref_level0.RData")
save(MM.FMT.level1.decon, file = "./Granulator_decon_MM_FMT_cycling_ref_level1.RData")

##### Level 0 annotation #####

nnls.res <- MM.FMT.decon$proportions$nnls_sig1
nnls.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

ols.res <- MM.FMT.decon$proportions$ols_sig1
ols.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

qprog.res <- MM.FMT.decon$proportions$qprog_sig1
qprog.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

qprogwc.res <- MM.FMT.decon$proportions$qprogwc_sig1
qprogwc.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

rls.res <- MM.FMT.decon$proportions$rls_sig1
rls.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

dtangle.res <- MM.FMT.decon$proportions$dtangle_sig1
dtangle.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

svr.res <- MM.FMT.decon$proportions$svr_sig1
svr.res$FMT <- factor(c(rep("FMT_OF", 7), rep("FMT_YF", 7)), levels = c("FMT_YF","FMT_OF"))

# GSEA analysis indicated differences in immune function 
wilcox.test(nnls.res$Immune[nnls.res$FMT == "FMT_YF"] , nnls.res$Immune[nnls.res$FMT == "FMT_OF"]) # p-value = 0.2639
wilcox.test(qprogwc.res$Immune[qprogwc.res$FMT == "FMT_YF"] , qprogwc.res$Immune[qprogwc.res$FMT == "FMT_OF"]) # p-value = 0.2639
wilcox.test(dtangle.res$Immune[dtangle.res$FMT == "FMT_YF"] , dtangle.res$Immune[dtangle.res$FMT == "FMT_OF"]) # p-value = 1

pdf(paste0(Sys.Date(),"_MM_FMT_ovary_immune_cells_NNLS.pdf"), height = 4, width = 5)
beeswarm(Immune ~ FMT, data = nnls.res, pch = 16, ylim = c(0,60), main = "NNLS deconvolution", pwcol = c(rep("#B1746FFF",7),rep("#800000FF",7)) )
text(1.5, 55, "ns")
dev.off()

pdf(paste0(Sys.Date(),"_MM_FMT_ovary_immune_cells_QPROGWC.pdf"), height = 4, width = 5)
beeswarm(Immune ~ FMT, data = qprogwc.res, pch = 16, ylim = c(0,20), main = "QPROGWC deconvolution", pwcol = c(rep("#B1746FFF",7),rep("#800000FF",7)) )
text(1.5, 15, "ns")
dev.off()

pdf(paste0(Sys.Date(),"_MM_FMT_ovary_immune_cells_DTANGLE.pdf"), height = 4, width = 5)
beeswarm(Immune ~ FMT, data = dtangle.res, pch = 16, ylim = c(0,10), main = "DTANGLE deconvolution", pwcol = c(rep("#B1746FFF",7),rep("#800000FF",7)) )
text(1.5, 8, "ns")
dev.off()

################################################################################
sink(file = paste(Sys.Date(),"_MM_FMT_ovarian_bulk_RNAseq_deconvolution_CSCDRNA_and_Granulator_session_Info.txt", sep =""))
sessionInfo()
sink()
