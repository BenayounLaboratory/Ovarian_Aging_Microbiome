options(stringsAsFactors = F)

library(qiime2R)
library(tidyverse)
library(ggplot2)
theme_set(theme_bw())  
library(ALDEx2)
library(dplyr)

###################################
# Menopause-microbiome project
# Functional abundance prediction assay - picrust2
# Aging cohort - males
###################################

###################################
# 1. Import picrust2 output files
###################################

# AC cohort - males

AC.males.16.pathway.pred <- read_qza("./q2-picrust2_output_AC16_males/pathway_abundance.qza")$data
AC.males.17.pathway.pred <- read_qza("./q2-picrust2_output_AC17_males/pathway_abundance.qza")$data
AC.males.22.pathway.pred <- read_qza("./q2-picrust2_output_AC22_males/pathway_abundance.qza")$data

# map file

map.file <- read.table("./metacyc_pathways_info.txt", sep = "\t")

# Convert values into integers

AC.males.16.pathway.pred <- round(AC.males.16.pathway.pred)
AC.males.17.pathway.pred <- round(AC.males.17.pathway.pred)
AC.males.22.pathway.pred <- round(AC.males.22.pathway.pred)

###################################
# 2. Perform ALDEx2 corr
###################################

set.seed(12345)

# AC cohort - males

# AC16

conds <- ifelse(grepl("4m", colnames(AC.males.16.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC16 = aldex(AC.males.16.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC16, 10)

# AC17

conds <- ifelse(grepl("4m", colnames(AC.males.17.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC17 = aldex(AC.males.17.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC17, 10)

# AC22

conds <- ifelse(grepl("4m", colnames(AC.males.22.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC22 = aldex(AC.males.22.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC22, 10)

# Calculate and correct p-values

# Identifying common Feature.IDs across all tables

common_functions <- Reduce(intersect, list(rownames(aldex2.AC16), rownames(aldex2.AC17), rownames(aldex2.AC22)))

length(common_functions)     #374

# Subsetting each table to only include the common features and then sorting
my.results.16.filt <- aldex2.AC16[rownames(aldex2.AC16) %in% common_functions, ]
my.results.17.filt <- aldex2.AC17[rownames(aldex2.AC17) %in% common_functions, ][order(rownames(my.results.16.filt)), ]
my.results.22.filt <- aldex2.AC17[rownames(aldex2.AC22) %in% common_functions, ][order(rownames(my.results.16.filt)), ]

# Filter based on directionality

combined.results.1 <- cbind(my.results.16.filt[4], my.results.17.filt[4], my.results.22.filt[4])
colnames(combined.results.1) <- c("diff.btw.16", "diff.btw.17", "diff.btw.22")

# Calculate the sign for each of the last four columns
signs <- sign(combined.results.1)

# Sum the signs for each row
row_sums <- rowSums(signs)

# Identify rows where the absolute sum is 4 (all signs are the same)
all_negative_rows <- rowSums(combined.results.1 < 0) == 3
all_positive_rows <- rowSums(combined.results.1 > 0) == 3

# Extract Feature.IDs where the condition is true
combined.results.negative <- combined.results.1[all_negative_rows, ]
combined.results.positive <- combined.results.1[all_positive_rows, ]

# Merge results & perform Fisher's method

combined.results.2 <- cbind(my.results.16.filt[c(10)], my.results.17.filt[10], my.results.22.filt[10])
colnames(combined.results.2) <- c("wi.ep.16", "wi.ep.17", "wi.ep.22")

results.negative <- combined.results.2[rownames(combined.results.2) %in% rownames(combined.results.negative), ]
results.positive <- combined.results.2[rownames(combined.results.2) %in% rownames(combined.results.positive), ]

# Create a vector of p-values from the four columns
p_values.neg <- results.negative[, c("wi.ep.16", "wi.ep.17", "wi.ep.22")]
p_values.pos <- results.positive[, c("wi.ep.16", "wi.ep.17", "wi.ep.22")]

# Multiply the p-values for each row
multiplied_p_values.neg <- apply(p_values.neg, 1, prod)
multiplied_p_values.pos <- apply(p_values.pos, 1, prod)

# Apply the Benjamini-Hochberg (BH) correction
bh_corrected_p_values.neg <- p.adjust(multiplied_p_values.neg, method = "BH")
bh_corrected_p_values.pos <- p.adjust(multiplied_p_values.pos, method = "BH")

# Create a new data frame with the BH-corrected p-values
result_df.neg <- data.frame(BH_Corrected_P_Value = bh_corrected_p_values.neg)
result_df.pos <- data.frame(BH_Corrected_P_Value = bh_corrected_p_values.pos)

# Print the result
print(result_df.neg)
print(result_df.pos)

# Count how many rows have BH-corrected p-values less than 0.05, 0.01, 0.005
count_significant.neg <- sum(result_df.neg$BH_Corrected_P_Value < 0.05)     # 25
count_significant.pos <- sum(result_df.pos$BH_Corrected_P_Value < 0.05)     # 34

results_df.AC.males <- rbind(result_df.neg, result_df.pos)

##############################################
# 3. Differential abundance analysis - combine effect size
##############################################

# Extracting the effect sizes for each feature from each cohort
effects <- data.frame(
  cohort16 = my.results.16.filt$effect,
  cohort17 = my.results.17.filt$effect,
  cohort22 = my.results.22.filt$effect)

rownames(effects) <- rownames(my.results.16.filt)

# Calculating SEM for each feature across the cohorts
effects$sem <- apply(effects[, -1], 1, function(x) sd(x) / sqrt(length(x)))

# The resulting 'effects' dataframe now contains the SEM for each feature in the 'sem' column
print(effects)

# Calculating the average effect size for each feature
effects$average_effect <- rowMeans(effects[,1:3])

# Add average effect size information to neg and pos tables above
# Convert row names to a column for both data frames
results_df.AC.males$Pathway <- rownames(results_df.AC.males)
effects$Pathway <- rownames(effects)

merged_df <- merge(results_df.AC.males, effects[,c("Pathway", "average_effect")], by = "Pathway", all.x = TRUE)

# Perform a left join to map descriptions based on the Pathway column.
merged_df <- merged_df %>% 
  left_join(map.file, by = c("Pathway" = "V1"))

# Optionally, rename the V2 column to "description".
colnames(merged_df)[colnames(merged_df) == "V2"] <- "description"

aldex2_filtered <- subset(merged_df, average_effect >= 1)
aldex2_filtered <- aldex2_filtered[!is.na(aldex2_filtered$description), ]

aldex2_filtered$adjP <- -log10(aldex2_filtered$BH_Corrected_P_Value)

aldex2_filtered_sorted <- aldex2_filtered[order(-aldex2_filtered$BH_Corrected_P_Value), ]

pdf(paste0(Sys.Date(), "_MM_AC_males_picrust2_barplot_down_in_OM.pdf"), width = 15, height = 5)
ggplot(aldex2_filtered_sorted, aes(x = reorder(description, adjP), y = adjP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + # Flip coordinates for horizontal bars
  labs(x = "description", y = "-log10(adjP)", title = "Pathway vs. Adjusted P-value") +
  theme_minimal()
dev.off()

aldex2_filtered <- subset(merged_df, average_effect <= -1)
aldex2_filtered <- aldex2_filtered[!is.na(aldex2_filtered$description), ]

aldex2_filtered$adjP <- -log10(aldex2_filtered$BH_Corrected_P_Value)

aldex2_filtered_sorted <- aldex2_filtered[order(-aldex2_filtered$adjP), ]

pdf(paste0(Sys.Date(), "_MM_AC_males_picrust2_barplot_up_in_OM.pdf"), width = 15, height = 5)
ggplot(aldex2_filtered_sorted, aes(x = reorder(description, adjP), y = adjP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + # Flip coordinates for horizontal bars
  labs(x = "description", y = "-log10(adjP)", title = "Pathway vs. Adjusted P-value") +
  theme_minimal()
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_AC_males_differential_abundance_analysis.txt", sep =""))
sessionInfo()
sink()


