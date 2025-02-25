setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Analysis_V7")

# Load necessary libraries
library(tidyverse)
library(limma)

###################################
# Menapoause-microbiome project
# Serum metabolomics data
# Perform PCA & generate plots
###################################

###################################
# 1. Load Data
###################################

################
# Load matrix & ID mapping data

data <- read.csv('/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Raw_data/fc_log_fdr_peak_ms2_filtered_matrix.csv', row.names = 1)
ID.mapping <- read.csv('/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Raw_data/ID_to_chemical_name_mapping.csv', row.names = 1)

################
# Load metadata and combine with matrix

metadata <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Raw_data/MM_metabolomics_serum_metadata.txt", header = TRUE, sep = "\t")

# Merge metadata with metabolomics data
data$Sample_ID <- rownames(data)
merged_data <- left_join(metadata, data, by = "Sample_ID")

# Convert Group and Batch to factors
merged_data$Group <- as.factor(merged_data$Group)
merged_data$Batch <- as.factor(merged_data$Batch)

# Split data by batch

data_batch1 <- merged_data[merged_data$Batch == 1, ]
data_batch2 <- merged_data[merged_data$Batch == 2, ]

# Transform data

t.data_batch1 <- as.data.frame(t(data_batch1))
t.data_batch2 <- as.data.frame(t(data_batch2))

###################################
# 2. Perform VSN normalization & run PCA
###################################

# VSN normalization
my.data.batch1.vsn <- as.data.frame(normalizeVSN(as.matrix(t.data_batch1[-c(1:3),])))
my.data.batch2.vsn <- as.data.frame(normalizeVSN(as.matrix(t.data_batch2[-c(1:3),])))

# PCA

pca.batch1 <- prcomp(t(my.data.batch1.vsn), scale. = TRUE)
pca.batch2 <- prcomp(t(my.data.batch2.vsn), scale. = TRUE)

pca.batch1.summary <- summary(pca.batch1)
pca.batch2.summary <- summary(pca.batch2)

# Subset metadata by batch
metadata_batch1 <- subset(metadata, Batch == 1)
metadata_batch2 <- subset(metadata, Batch == 2)

# Create a color vector: purple for FMT_YF, orange for FMT_EF
colors_batch1 <- ifelse(metadata_batch1$Group == "FMT_YF", "purple", "orange")
colors_batch2 <- ifelse(metadata_batch2$Group == "FMT_YF", "purple", "orange")

# Save PCA plot
pdf(paste0(Sys.Date(), "_MM_Metabolomics_batch1_PCA_plot.pdf"))
plot(pca.batch1.summary$x[, 1], pca.batch1.summary$x[, 2],
     col = colors_batch1,
     pch = 16,
     cex = 2,
     xlab = paste("PC1 (", round(summary(pca.result)$importance[2, 1] * 100, 1), "%)", sep = ""),
     ylab = paste("PC2 (", round(summary(pca.result)$importance[2, 2] * 100, 1), "%)", sep = ""),
     main = "Batch 1 PCA: FMT_YF vs. FMT_EF")
legend("topright", legend = names(group_colors), col = group_colors, pch = 16)
dev.off()

pdf(paste0(Sys.Date(), "_MM_Metabolomics_batch2_PCA_plot.pdf"))
plot(pca.batch2.summary$x[, 1], pca.batch2.summary$x[, 2],
     col = colors_batch2,
     pch = 16,
     cex = 2,
     xlab = paste("PC1 (", round(summary(pca.result)$importance[2, 1] * 100, 1), "%)", sep = ""),
     ylab = paste("PC2 (", round(summary(pca.result)$importance[2, 2] * 100, 1), "%)", sep = ""),
     main = "Batch 2 PCA: FMT_YF vs. FMT_EF")
legend("topright", legend = names(group_colors), col = group_colors, pch = 16)
dev.off()

################################################################################
# Save session info
sink(file = paste0(Sys.Date(), "_MM_metabolomics_analysis_PCA_Session_Info.txt"))
sessionInfo()
sink()
