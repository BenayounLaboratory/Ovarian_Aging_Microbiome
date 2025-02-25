setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Reviewer_panel_metagenomics_hormone_correlation")

library(ggcorrplot)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(dplyr)

#############################################################
# Menopause-microbiome project - revision
# Metagenomics - hormone data correlation analysis
# Use CLR-transformed dataset
#############################################################

#############################################################
# 1. Load data
#############################################################

# FMT metagenomics data & metadata

metagenomics.clr <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2025-01-21/2025-01-30_MM_FMT_metagenomics_CLR_transformed_data.txt", header = TRUE, sep = "\t")

sig.species <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2025-01-21/2025-01-30_MM_FMT_metagenomics_sig_features.txt", header = TRUE, sep = "\t")
colnames(sig.species) <- c("Species")

FMT.metadata <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2024-12-02/FMT_Metagenomics_metadata.txt", 
                           header = TRUE)

# Hormone data

my.hormone.data <- read.table("./MM_FMT_hormone_data_updated_20250117.txt", header = TRUE)

# Update hormone data - maintain overlapping data with metagenomics only

my.hormone.data.c1 <- my.hormone.data[my.hormone.data$Cohort == "FMT_1",]
my.hormone.data.c1$Metagenomics_ID <- FMT.metadata$Sample_ID

#############################################################
# 2. Filter metagenomics data - significant species only (top & bottom 10 species)
#############################################################

# Extract significant species

metagenomics.sig.data <- metagenomics.clr[sig.species$Species,]

# Combine with hormone data

metagenomics.sig.data.transposed <- as.data.frame(t(metagenomics.sig.data))

metagenomics.sig.data.transposed$Metagenomics_ID <- rownames(metagenomics.sig.data.transposed)

combined_data <- merge(my.hormone.data.c1, metagenomics.sig.data.transposed, by="Metagenomics_ID")

write.table(combined_data, file = paste0(Sys.Date(), "_MM_revision_FMT_cohort_metagenomics_sig_20_species_and_hormone_data_combined.txt"), quote = FALSE, sep = "\t")

#############################################################
# 3. Compute Spearman correlation and p-values
#############################################################

# Define species and hormone columns

hormone_cols <- c("AMH", "FSH_Multiplex", "INHBA")
species_cols <- colnames(combined_data)[!(colnames(combined_data) %in% c("Sample_ID", "Mouse_ID", "Treatment", "Cohort", "Metagenomics_ID", "FSH_US", hormone_cols))]

# Convert non-numeric columns (hormones & species) to numeric

combined_data[, c(species_cols, hormone_cols)] <- lapply(combined_data[, c(species_cols, hormone_cols)], function(x) {
  if (is.character(x) | is.factor(x)) as.numeric(as.character(x)) else x
})

# Initialize - store correlation coefficients and p-values

cor_matrix <- matrix(NA, nrow=length(species_cols), ncol=length(hormone_cols), dimnames=list(species_cols, hormone_cols))
pval_matrix <- matrix(NA, nrow=length(species_cols), ncol=length(hormone_cols), dimnames=list(species_cols, hormone_cols))

# Compute Spearman correlation and p-values

for (microbe in species_cols) {
  for (hormone in hormone_cols) {
    test_result <- cor.test(combined_data[[microbe]], combined_data[[hormone]], method = "spearman")
    cor_matrix[microbe, hormone] <- test_result$estimate 
    pval_matrix[microbe, hormone] <- test_result$p.value 
  }
}

# Format p-values for display

pval_labels <- apply(pval_matrix, c(1,2), function(x) formatC(x, format = "e", digits = 2)) 

# Generate heatmap
pdf(paste0(Sys.Date(), "_MM_revision_FMT_cohort_metagenomics_sig_20_species_and_hormone_data_correlation_heatmap.pdf"), width = 10, height = 8)
pheatmap(cor_matrix, 
         cluster_rows = TRUE, cluster_cols = TRUE,  
         display_numbers = pval_labels, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation: Species & Hormones",
         fontsize = 10)
dev.off()

# Save correlation and p-value tables

write.table(cor_matrix, file = paste0(Sys.Date(), "_MM_revision_FMT_combined_correlation_coefficients.txt"), sep = "\t", quote = FALSE)
write.table(pval_matrix, file = paste0(Sys.Date(), "_MM_revision_FMT_combined_pvalues.txt"), sep = "\t", quote = FALSE)

#############################################################
# 4. By-group analysis
#############################################################

# Separate data by treatment
combined_YF <- combined_data[combined_data$Treatment == "FMT_YF", ]
combined_EF <- combined_data[combined_data$Treatment == "FMT_EF", ]

##########################
# Compute Correlation for FMT-YF
##########################

# Initialize

cor_matrix_YF <- matrix(NA, nrow=length(species_cols), ncol=length(hormone_cols), dimnames=list(species_cols, hormone_cols))
pval_matrix_YF <- matrix(NA, nrow=length(species_cols), ncol=length(hormone_cols), dimnames=list(species_cols, hormone_cols))

# Compute correlations for FMT-YF

for (microbe in species_cols) {
  for (hormone in hormone_cols) {
    test_result <- cor.test(combined_YF[[microbe]], combined_YF[[hormone]], method = "spearman")
    cor_matrix_YF[microbe, hormone] <- test_result$estimate
    pval_matrix_YF[microbe, hormone] <- test_result$p.value
  }
}

##########################
# Compute Correlation for FMT-EF
##########################

# Initialize matrices

cor_matrix_EF <- matrix(NA, nrow=length(species_cols), ncol=length(hormone_cols), dimnames=list(species_cols, hormone_cols))
pval_matrix_EF <- matrix(NA, nrow=length(species_cols), ncol=length(hormone_cols), dimnames=list(species_cols, hormone_cols))

# Compute correlations for FMT-EF

for (microbe in species_cols) {
  for (hormone in hormone_cols) {
    test_result <- cor.test(combined_EF[[microbe]], combined_EF[[hormone]], method = "spearman")
    cor_matrix_EF[microbe, hormone] <- test_result$estimate
    pval_matrix_EF[microbe, hormone] <- test_result$p.value
  }
}

##########################
# Generate heatmaps
##########################

pdf(paste0(Sys.Date(), "_MM_revision_FMT_cohort_metagenomics_sig_20_species_and_hormone_data_correlation_heatmap_FMT-YF_only.pdf"), width = 10, height = 8)
pheatmap(cor_matrix_YF, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = apply(pval_matrix_YF, c(1,2), function(x) formatC(x, format = "e", digits = 2)),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation (FMT-YF Only)",
         fontsize = 10)
dev.off()

pdf(paste0(Sys.Date(), "_MM_revision_FMT_cohort_metagenomics_sig_20_species_and_hormone_data_correlation_heatmap_FMT-EF_only.pdf"), width = 10, height = 8)
pheatmap(cor_matrix_EF, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = apply(pval_matrix_EF, c(1,2), function(x) formatC(x, format = "e", digits = 2)),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation (FMT-EF Only)",
         fontsize = 10)
dev.off()

##########################
# Save tables
##########################

write.table(cor_matrix_YF, file = paste0(Sys.Date(), "_MM_revision_FMT_YF_correlation_coefficients.txt"), sep = "\t", quote = FALSE)
write.table(pval_matrix_YF, file = paste0(Sys.Date(), "_MM_revision_FMT_YF_pvalues.txt"), sep = "\t", quote = FALSE)

write.table(cor_matrix_EF, file = paste0(Sys.Date(), "_MM_revision_FMT_EF_correlation_coefficients.txt"), sep = "\t", quote = FALSE)
write.table(pval_matrix_EF, file = paste0(Sys.Date(), "_MM_revision_FMT_EF_pvalues.txt"), sep = "\t", quote = FALSE)

#############################################################
sink(file = paste0(Sys.Date(),"_MM_revision_FMT_cohort_metagenomics_hormone_correlation_analysis_session_Info.txt"))
sessionInfo()
sink()
