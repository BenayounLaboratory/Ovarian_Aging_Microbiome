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
# AC cohort - females
###################################

##############################################
# 1. Import data files
##############################################

# metadata

my.results.16  <- read_qza("./aldex2_AC_females_cohort16/differentials.qza")$data
my.results.17  <- read_qza("./aldex2_AC_females_cohort17/differentials.qza")$data
my.results.22  <- read_qza("./aldex2_AC_females_cohort22/differentials.qza")$data
my.results.24  <- read_qza("./aldex2_AC_females_cohort24/differentials.qza")$data

my.taxonomy.16 <- read_qza("../../5_Taxonomy/silva_MM_AC16_females_taxonomy.qza")$data
my.taxonomy.17 <- read_qza("../../5_Taxonomy/silva_MM_AC17_females_taxonomy.qza")$data
my.taxonomy.22 <- read_qza("../../5_Taxonomy/silva_MM_AC22_females_taxonomy.qza")$data
my.taxonomy.24 <- read_qza("../../5_Taxonomy/silva_MM_AC24_females_taxonomy.qza")$data

my.taxonomy.combined <- rbind(my.taxonomy.16, my.taxonomy.17, my.taxonomy.22, my.taxonomy.24)

my.taxonomy.combined.cl <- my.taxonomy.combined %>%
  distinct(Feature.ID, .keep_all = TRUE)

##############################################
# 2. Differential abundance analysis - calculate and correct p-values
##############################################

# Identifying common Feature.IDs across all tables
common_features <- Reduce(intersect, list(my.results.16$Feature.ID, my.results.17$Feature.ID, 
                                          my.results.22$Feature.ID, my.results.24$Feature.ID))

length(common_features)     #299

# Subsetting each table to only include the common features and then sorting
my.results.16.filt <- my.results.16[my.results.16$Feature.ID %in% common_features, ][order(my.results.16$Feature.ID[my.results.16$Feature.ID %in% common_features]), ]
my.results.17.filt <- my.results.17[my.results.17$Feature.ID %in% common_features, ][order(my.results.17$Feature.ID[my.results.17$Feature.ID %in% common_features]), ]
my.results.22.filt <- my.results.22[my.results.22$Feature.ID %in% common_features, ][order(my.results.22$Feature.ID[my.results.22$Feature.ID %in% common_features]), ]
my.results.24.filt <- my.results.24[my.results.24$Feature.ID %in% common_features, ][order(my.results.24$Feature.ID[my.results.24$Feature.ID %in% common_features]), ]

# Filter based on directionality

combined.results.1 <- cbind(my.results.16.filt[c(1,5)], my.results.17.filt[5], my.results.22.filt[5], my.results.24.filt[5])
colnames(combined.results.1) <- c("Feature.ID", "diff.btw.16", "diff.btw.17", "diff.btw.22", "diff.btw.24")

# Calculate the sign for each of the last four columns
signs <- sign(combined.results.1[, 2:5])

# Sum the signs for each row
row_sums <- rowSums(signs)

# Identify rows where the absolute sum is 4 (all signs are the same)
all_negative_rows <- rowSums(combined.results.1[, 2:5] < 0) == 4
all_positive_rows <- rowSums(combined.results.1[, 2:5] > 0) == 4

# Extract Feature.IDs where the condition is true
combined.results.negative <- combined.results.1[all_negative_rows, ]
combined.results.positive <- combined.results.1[all_positive_rows, ]

# Merge results & perform Fisher's method

combined.results.2 <- cbind(my.results.16.filt[c(1,11)], my.results.17.filt[11], my.results.22.filt[11], my.results.24.filt[11])
colnames(combined.results.2) <- c("Feature.ID", "wi.ep.16", "wi.ep.17", "wi.ep.22", "wi.ep.24")

results.negative <- combined.results.2[combined.results.2$Feature.ID %in% combined.results.negative$Feature.ID, ]
results.positive <- combined.results.2[combined.results.2$Feature.ID %in% combined.results.positive$Feature.ID, ]

# Create a vector of p-values from the four columns
p_values.neg <- results.negative[, c("wi.ep.16", "wi.ep.17", "wi.ep.22", "wi.ep.24")]
p_values.pos <- results.positive[, c("wi.ep.16", "wi.ep.17", "wi.ep.22", "wi.ep.24")]

# Multiply the p-values for each row
multiplied_p_values.neg <- apply(p_values.neg, 1, prod)
multiplied_p_values.pos <- apply(p_values.pos, 1, prod)

# Apply the Benjamini-Hochberg (BH) correction
bh_corrected_p_values.neg <- p.adjust(multiplied_p_values.neg, method = "BH")
bh_corrected_p_values.pos <- p.adjust(multiplied_p_values.pos, method = "BH")

# Create a new data frame with the BH-corrected p-values
result_df.neg <- data.frame(Feature.ID = results.negative$Feature.ID, BH_Corrected_P_Value = bh_corrected_p_values.neg)
result_df.pos <- data.frame(Feature.ID = results.positive$Feature.ID, BH_Corrected_P_Value = bh_corrected_p_values.pos)

# Print the result
print(result_df.neg)
print(result_df.pos)

write.table(result_df.neg, file = "BH_Corrected_p_values_DOWN_with_age.txt", quote = FALSE, row.names = FALSE)
write.table(result_df.pos, file = "BH_Corrected_p_values_UP_with_age.txt", quote = FALSE, row.names = FALSE)

# Count how many rows have BH-corrected p-values less than 0.05, 0.01, 0.005
count_significant.neg <- sum(result_df.neg$BH_Corrected_P_Value < 0.05)     # 41
count_significant.pos <- sum(result_df.pos$BH_Corrected_P_Value < 0.05)     # 23

results_df.neg.pval.0.05 <- result_df.neg[result_df.neg$BH_Corrected_P_Value < 0.05,]
results_df.pos.pval.0.05 <- result_df.pos[result_df.pos$BH_Corrected_P_Value < 0.05,]

write.table(results_df.neg.pval.0.05, file = "BH_Corrected_p_values_DOWN_with_age_pval_0.05.txt", quote = FALSE, row.names = FALSE)
write.table(results_df.pos.pval.0.05, file = "BH_Corrected_p_values_UP_with_age_pval_0.05.txt", quote = FALSE, row.names = FALSE)

results.combined.pval.0.05 <- rbind(results_df.neg.pval.0.05, results_df.pos.pval.0.05)

##############################################
# 3. Differential abundance analysis - combine effect size
##############################################

# Extracting the effect sizes for each feature from each cohort
effects <- data.frame(
  feature.ID = my.results.16.filt$Feature.ID, # Assuming the feature IDs are aligned
  cohort16 = my.results.16.filt$effect,
  cohort17 = my.results.17.filt$effect,
  cohort22 = my.results.22.filt$effect,
  cohort24 = my.results.24.filt$effect
)

# Calculating SEM for each feature across the cohorts
effects$sem <- apply(effects[, -1], 1, function(x) sd(x) / sqrt(length(x)))

# The resulting 'effects' dataframe now contains the SEM for each feature in the 'sem' column
print(effects)

# Calculating the average effect size for each feature
effects$average_effect <- rowMeans(effects[,2:5])

# Add average effect size information to neg and pos tables above
effects$Feature.ID <- effects$feature.ID

results.combined.pval.0.05 <- merge(results.combined.pval.0.05, effects[, c("Feature.ID", "average_effect")], by = "Feature.ID", all.x = TRUE)

# Filter siginificant features
results.combined.pval.0.05$signif <- with(results.combined.pval.0.05, average_effect > 1 | average_effect < -1)

significant_effects <- results.combined.pval.0.05[results.combined.pval.0.05$signif == TRUE,]
significant_effects <- significant_effects[order(significant_effects$average_effect < 0, significant_effects$BH_Corrected_P_Value), ]

# Plotting heatmap of significant features

significant_effects_unique <- significant_effects[!duplicated(significant_effects$Feature.ID) & !duplicated(significant_effects$Feature.ID, fromLast = TRUE), ]

significant_effects_unique <- significant_effects[!duplicated(significant_effects$Feature.ID) & !duplicated(significant_effects$Feature.ID, fromLast = TRUE), ]

significant_effects_unique$minus_log10_pvalue <- -log10(significant_effects_unique$BH_Corrected_P_Value)

# Rename Feature ID with taxon name

# Trim the Taxon information and format it by removing "D_#__" and handling empty values
my.taxonomy.combined.cl$Taxon_formatted <- sapply(strsplit(as.character(my.taxonomy.combined.cl$Taxon), ";"), function(taxon) {
  # Remove "D_#__" from each element and replace empty values with NA
  formatted_taxon <- gsub("D_\\d+__", "", taxon)
  formatted_taxon <- ifelse(formatted_taxon == "", "NA", formatted_taxon)
  # Keep only up to D_7 (first 8 elements)
  formatted_taxon <- formatted_taxon[1:min(length(formatted_taxon), 8)]
  # Collapse back into a single string, replacing empty slots with "NA"
  paste(formatted_taxon, collapse = ";")
})

# Create a lookup table with Feature.ID and the formatted Taxon
lookup_table <- data.frame(Feature.ID = my.taxonomy.combined.cl$Feature.ID, Taxon = my.taxonomy.combined.cl$Taxon_formatted)

# Replace Feature.ID in my.results.VCD with the corresponding Taxon from the lookup table
significant_effects_unique$Feature.ID <- sapply(significant_effects_unique$Feature.ID, function(id) {
  if(id %in% lookup_table$Feature.ID) {
    return(lookup_table$Taxon[lookup_table$Feature.ID == id])
  } else {
    return(id) # Return the original ID if not found in the lookup table
  }
})

# Creating the bubble plot
ggplot(significant_effects_unique, aes(x = reorder(Feature.ID, average_effect), y = average_effect, size = minus_log10_pvalue)) +
  geom_point(aes(color = average_effect)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 10), name = "-log10(p-value)") + # Use a meaningful size scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank()) +
  labs(title = "Significant Features - AC females", 
       x = "Features", 
       y = "Effect", 
       size = "-log10(P-value)", 
       color = "Effect")

ggsave("MM_AC_females_Differential_bubble_plot.pdf", height=15, width=7, device="pdf", useDingbats=F)


###########################################################
sink(file = paste(Sys.Date(),"MM_AC_females_differential_abundance_analysis.txt", sep =""))
sessionInfo()
sink()
