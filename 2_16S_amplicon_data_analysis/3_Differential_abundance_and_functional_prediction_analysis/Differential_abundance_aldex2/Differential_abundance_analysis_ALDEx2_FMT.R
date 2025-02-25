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
# FMT cohort
###################################

##############################################
# 1. Import data files
##############################################

# metadata

metadata.FMT <- read_q2metadata("../../1_Manifest_and_metadata/MM_FMT_AC_PostFMT_sample-metadata.tsv")
metadata.FMT$fmt <- factor(FMT.metadata$fmt, levels = c("FMT_Y", "FMT_O"))

SVs.FMT      <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_FMT_AC_PostFMT_table.qza")$data
results.FMT  <- read_qza("../../6_Differential_abundance_and_functional_prediction_analysis/Differential_abundance_aldex2/aldex2_FMT/differentials.qza")$data

taxonomy.FMT <- read_qza("../../5_Taxonomy/silva_MM_FMT_AC_PostFMT_taxonomy.qza")$data
tree.FMT     <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_VCD_rooted-tree.qza")$data

##############################################
# 2. Differential abundance analysis
##############################################

# Rename Feature ID with taxon name

# Trim the Taxon information and format it by removing "D_#__" and handling empty values
taxonomy.FMT$Taxon_formatted <- sapply(strsplit(as.character(taxonomy.FMT$Taxon), ";"), function(taxon) {
  # Remove "D_#__" from each element and replace empty values with NA
  formatted_taxon <- gsub("D_\\d+__", "", taxon)
  formatted_taxon <- ifelse(formatted_taxon == "", "NA", formatted_taxon)
  # Keep only up to D_7 (first 8 elements)
  formatted_taxon <- formatted_taxon[1:min(length(formatted_taxon), 8)]
  # Collapse back into a single string, replacing empty slots with "NA"
  paste(formatted_taxon, collapse = ";")
})

# Create a lookup table with Feature.ID and the formatted Taxon
lookup_table <- data.frame(Feature.ID = taxonomy.FMT$Feature.ID, Taxon = taxonomy.FMT$Taxon_formatted)

# Replace Feature.ID in results.FMT with the corresponding Taxon from the lookup table
results.FMT$Feature.ID <- sapply(results.FMT$Feature.ID, function(id) {
  if(id %in% lookup_table$Feature.ID) {
    return(lookup_table$Taxon[lookup_table$Feature.ID == id])
  } else {
    return(id) # Return the original ID if not found in the lookup table
  }
})

# Divide neg vs. pos diff.btw entries

result_df.neg <- subset(results.FMT, diff.btw < 0)
result_df.pos <- subset(results.FMT, diff.btw > 0)

# Print the result
print(result_df.neg)
print(result_df.pos)

write.table(result_df.neg, file = "DOWN_with_FMT.txt", quote = FALSE, row.names = FALSE)
write.table(result_df.pos, file = "UP_with_FMT.txt", quote = FALSE, row.names = FALSE)

# Count how many rows have BH-corrected p-values less than 0.05, 0.01, 0.005
count_significant.neg <- sum(result_df.neg$wi.ep < 0.05)     # 19
count_significant.pos <- sum(result_df.pos$wi.ep < 0.05)     # 24

results_df.neg.pval.0.05 <- result_df.neg[result_df.neg$wi.ep < 0.05,]
results_df.pos.pval.0.05 <- result_df.pos[result_df.pos$wi.ep < 0.05,]

write.table(results_df.neg.pval.0.05, file = "DOWN_with_FMT_pval_0.05.txt", quote = FALSE, row.names = FALSE)
write.table(results_df.pos.pval.0.05, file = "UP_with_FMT_pval_0.05.txt", quote = FALSE, row.names = FALSE)

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
  labs(title = "Significant Features - FMT", 
       x = "Features", 
       y = "Effect", 
       size = "-log10(P-value)", 
       color = "Effect")

ggsave("MM_FMT_Differential_bubble_plot.pdf", height=15, width=7, device="pdf", useDingbats=F)

##############################################
# 6. Plot single features
##############################################

clr <- apply(log2(SVs.FMT+0.5), 2, function(x) x-mean(x))
clr %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key=SampleID, value=CLR) %>%
  filter(Feature.ID=="d07ab1972ef6ebba6f2f95d9500608ae") %>%
  left_join(metadata.FMT) %>%
  ggplot(aes(x=fmt, y=CLR, fill=fmt)) +
  stat_summary(geom="bar", color="black") +
  geom_jitter(width=0.2, height=0, shape=21) +
  theme_q2r() +
  theme(legend.position="none")
ggsave("d07ab1972ef6ebba6f2f95d9500608ae_aldexbar.pdf", height=2, width=1.5, device="pdf") 

clr %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key=SampleID, value=CLR) %>%
  filter(Feature.ID=="793929587b93e71a3b11e9cd667bea5e") %>%
  left_join(metadata.FMT) %>%
  ggplot(aes(x=fmt, y=CLR, fill=fmt)) +
  stat_summary(geom="bar", color="black") +
  geom_jitter(width=0.2, height=0, shape=21) +
  theme_q2r() +
  theme(legend.position="none")
ggsave("793929587b93e71a3b11e9cd667bea5e", height=2, width=1.5, device="pdf") 


###########################################################
sink(file = paste(Sys.Date(),"MM_FMT_differential_abundance_analysis.txt", sep =""))
sessionInfo()
sink()
