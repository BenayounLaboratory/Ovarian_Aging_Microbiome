library('qiime2R')
library(tidyverse)
library('beeswarm')
library(ggrepel) # for offset labels
library(ggtree) # for visualizing phylogenetic trees
library(ape) # for manipulating phylogenetic trees
library('ggplot2')

###################################
# Menopause-microbiome project
# Differenial abundance analysis - aldex2
# VCD cohort
###################################

##############################################
# 1. Import data files
##############################################

# metadata

metadata.VCD <- read_q2metadata("../../1_Manifest_and_metadata/MM_VCD_cohort_sample-metadata.tsv")
metadata.VCD$VCD <- factor(VCD.metadata$treatment, levels = c("CTL", "VCD"))

SVs.VCD      <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_VCD_table.qza")$data
results.VCD  <- read_qza("../../6_Differential_abundance_and_functional_prediction_analysis/Differential_abundance_aldex2/aldex2_VCD/differentials.qza")$data

taxonomy.VCD <- read_qza("../../5_Taxonomy/silva_MM_VCD_taxonomy.qza")$data
tree.VCD     <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_VCD_rooted-tree.qza")$data

##############################################
# 2. Differential abundance analysis
##############################################

# Rename Feature ID with taxon name

# Trim the Taxon information and format it by removing "D_#__" and handling empty values
taxonomy.VCD$Taxon_formatted <- sapply(strsplit(as.character(taxonomy.VCD$Taxon), ";"), function(taxon) {
  # Remove "D_#__" from each element and replace empty values with NA
  formatted_taxon <- gsub("D_\\d+__", "", taxon)
  formatted_taxon <- ifelse(formatted_taxon == "", "NA", formatted_taxon)
  # Keep only up to D_7 (first 8 elements)
  formatted_taxon <- formatted_taxon[1:min(length(formatted_taxon), 8)]
  # Collapse back into a single string, replacing empty slots with "NA"
  paste(formatted_taxon, collapse = ";")
})

# Create a lookup table with Feature.ID and the formatted Taxon
lookup_table <- data.frame(Feature.ID = taxonomy.VCD$Feature.ID, Taxon = taxonomy.VCD$Taxon_formatted)

# Replace Feature.ID in results.VCD with the corresponding Taxon from the lookup table
results.VCD$Feature.ID <- sapply(results.VCD$Feature.ID, function(id) {
  if(id %in% lookup_table$Feature.ID) {
    return(lookup_table$Taxon[lookup_table$Feature.ID == id])
  } else {
    return(id) # Return the original ID if not found in the lookup table
  }
})

# Divide neg vs. pos diff.btw entries

result_df.neg <- subset(results.VCD, diff.btw < 0)
result_df.pos <- subset(results.VCD, diff.btw > 0)

# Print the result
print(result_df.neg)
print(result_df.pos)

write.table(result_df.neg, file = "DOWN_with_VCD.txt", quote = FALSE, row.names = FALSE)
write.table(result_df.pos, file = "UP_with_VCD.txt", quote = FALSE, row.names = FALSE)

# Count how many rows have BH-corrected p-values less than 0.05, 0.01, 0.005
count_significant.neg <- sum(result_df.neg$wi.ep < 0.05)     # 9
count_significant.pos <- sum(result_df.pos$wi.ep < 0.05)     # 20

results_df.neg.pval.0.05 <- result_df.neg[result_df.neg$wi.ep < 0.05,]
results_df.pos.pval.0.05 <- result_df.pos[result_df.pos$wi.ep < 0.05,]

write.table(results_df.neg.pval.0.05, file = "DOWN_with_VCD_pval_0.05.txt", quote = FALSE, row.names = FALSE)
write.table(results_df.pos.pval.0.05, file = "UP_with_VCD_pval_0.05.txt", quote = FALSE, row.names = FALSE)

results.combined.pval.0.05 <- rbind(results_df.neg.pval.0.05, results_df.pos.pval.0.05)

results.combined.pval.0.05$signif <- with(results.combined.pval.0.05, effect > 1 | effect < -1)

# Filtering for significant features
significant_effects <- results.combined.pval.0.05[results.combined.pval.0.05$signif == TRUE,]
significant_effects <- significant_effects[order(significant_effects$effect < 0, significant_effects$wi.ep), ]

# Plotting heatmap of significant features
significant_effects_unique <- significant_effects[!duplicated(significant_effects$Feature.ID) & !duplicated(significant_effects$Feature.ID, fromLast = TRUE), ]

significant_effects_unique$minus_log10_pvalue <- -log10(significant_effects_unique$wi.ep)

# Creating the bubble plot
ggplot(significant_effects_unique, aes(x = reorder(Feature.ID, effect), y = effect, size = minus_log10_pvalue)) +
  geom_point(aes(color = effect)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 10), name = "-log10(p-value)") + # Use a meaningful size scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) +
  labs(title = "Significant Features - VCD", 
       x = "Features", 
       y = "Effect", 
       size = "-log10(P-value)", 
       color = "Effect")

ggsave("MM_VCD_Differential_bubble_plot.pdf", height=15, width=7, device="pdf", useDingbats=F)

###########################################################
sink(file = paste(Sys.Date(),"MM_VCD_differential_abundance_analysis.txt", sep =""))
sessionInfo()
sink()
