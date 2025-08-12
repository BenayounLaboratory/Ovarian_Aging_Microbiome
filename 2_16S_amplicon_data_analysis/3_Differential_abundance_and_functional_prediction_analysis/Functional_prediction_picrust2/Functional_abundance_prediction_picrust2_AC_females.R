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
# Aging cohort - females
###################################

###################################
# 1. Import picrust2 output files
###################################

# AC cohort - females

AC.females.16.pathway.pred <- read_qza("./q2-picrust2_output_AC16_females/pathway_abundance.qza")$data
AC.females.17.pathway.pred <- read_qza("./q2-picrust2_output_AC17_females/pathway_abundance.qza")$data
AC.females.22.pathway.pred <- read_qza("./q2-picrust2_output_AC22_females/pathway_abundance.qza")$data
AC.females.24.pathway.pred <- read_qza("./q2-picrust2_output_AC24_females/pathway_abundance.qza")$data

# map file

map.file <- read.table("./metacyc_pathways_info.txt", sep = "\t")

# Convert values into integers

AC.females.16.pathway.pred <- round(AC.females.16.pathway.pred)
AC.females.17.pathway.pred <- round(AC.females.17.pathway.pred)
AC.females.22.pathway.pred <- round(AC.females.22.pathway.pred)
AC.females.24.pathway.pred <- round(AC.females.24.pathway.pred)

###################################
# 2. Perform ALDEx2 corr
###################################

set.seed(12345)

# AC cohort - females

# AC16

conds <- ifelse(grepl("4m", colnames(AC.females.16.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC16 = aldex(AC.females.16.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC16, 10)

# AC17

conds <- ifelse(grepl("4m", colnames(AC.females.17.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC17 = aldex(AC.females.17.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC17, 10)

# AC22

conds <- ifelse(grepl("4m", colnames(AC.females.22.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC22 = aldex(AC.females.22.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC22, 10)

# AC24

conds <- ifelse(grepl("4m", colnames(AC.females.24.pathway.pred)), "4m", "20m")
conds.int <- ifelse(grepl("4m", conds), "1", "2")

# Run Aldex2

aldex2.AC24 = aldex(AC.females.24.pathway.pred, conds, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2.AC24, 10)

# Calculate and correct p-values

# Identifying common Feature.IDs across all tables

common_functions <- Reduce(intersect, list(rownames(aldex2.AC16), rownames(aldex2.AC17), rownames(aldex2.AC22), rownames(aldex2.AC24)))

length(common_functions)     #371

# Subsetting each table to only include the common features and then sorting
my.results.16.filt <- aldex2.AC16[rownames(aldex2.AC16) %in% common_functions, ]
my.results.17.filt <- aldex2.AC17[rownames(aldex2.AC17) %in% common_functions, ][order(rownames(my.results.16.filt)), ]
my.results.22.filt <- aldex2.AC17[rownames(aldex2.AC22) %in% common_functions, ][order(rownames(my.results.16.filt)), ]
my.results.24.filt <- aldex2.AC17[rownames(aldex2.AC24) %in% common_functions, ][order(rownames(my.results.16.filt)), ]

# Filter based on directionality

combined.results.1 <- cbind(my.results.16.filt[4], my.results.17.filt[4], my.results.22.filt[4], my.results.24.filt[4])
colnames(combined.results.1) <- c("diff.btw.16", "diff.btw.17", "diff.btw.22", "diff.btw.24")

# Calculate the sign for each of the last four columns
signs <- sign(combined.results.1)

# Sum the signs for each row
row_sums <- rowSums(signs)

# Identify rows where the absolute sum is 4 (all signs are the same)
all_negative_rows <- rowSums(combined.results.1 < 0) == 4
all_positive_rows <- rowSums(combined.results.1 > 0) == 4

# Extract Feature.IDs where the condition is true
combined.results.negative <- combined.results.1[all_negative_rows, ]
combined.results.positive <- combined.results.1[all_positive_rows, ]

# Merge results & perform Fisher's method

combined.results.2 <- cbind(my.results.16.filt[c(10)], my.results.17.filt[10], my.results.22.filt[10], my.results.24.filt[10])
colnames(combined.results.2) <- c("wi.ep.16", "wi.ep.17", "wi.ep.22", "wi.ep.24")

results.negative <- combined.results.2[rownames(combined.results.2) %in% rownames(combined.results.negative), ]
results.positive <- combined.results.2[rownames(combined.results.2) %in% rownames(combined.results.positive), ]

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
result_df.neg <- data.frame(BH_Corrected_P_Value = bh_corrected_p_values.neg)
result_df.pos <- data.frame(BH_Corrected_P_Value = bh_corrected_p_values.pos)

# Print the result
print(result_df.neg)
print(result_df.pos)

# Count how many rows have BH-corrected p-values less than 0.05, 0.01, 0.005
count_significant.neg <- sum(result_df.neg$BH_Corrected_P_Value < 0.05)     # 19
count_significant.pos <- sum(result_df.pos$BH_Corrected_P_Value < 0.05)     # 32

results_df.AC.females <- rbind(result_df.neg, result_df.pos)

##############################################
# 3. Differential abundance analysis - combine effect size
##############################################

# AC cohort - females

# Extracting the effect sizes for each feature from each cohort
effects <- data.frame(
  cohort16 = my.results.16.filt$effect,
  cohort17 = my.results.17.filt$effect,
  cohort22 = my.results.22.filt$effect,
  cohort24 = my.results.24.filt$effect
)

rownames(effects) <- rownames(my.results.16.filt)

# Calculating SEM for each feature across the cohorts
effects$sem <- apply(effects[, -1], 1, function(x) sd(x) / sqrt(length(x)))

# The resulting 'effects' dataframe now contains the SEM for each feature in the 'sem' column
print(effects)

# Calculating the average effect size for each feature
effects$average_effect <- rowMeans(effects[,1:4])

# Add average effect size information to neg and pos tables above
results_df.AC.females$Pathway <- rownames(results_df.AC.females)
effects$Pathway <- rownames(effects)

merged_df <- merge(results_df.AC.females, effects[,c("Pathway", "average_effect")], by = "Pathway", all.x = TRUE)

merged_df <- merged_df %>% 
  left_join(map.file, by = c("Pathway" = "V1"))

# Rename the V2 column to "description".
colnames(merged_df)[colnames(merged_df) == "V2"] <- "description"

aldex2_filtered <- subset(merged_df, average_effect >= 0.5)
aldex2_filtered <- aldex2_filtered[!is.na(aldex2_filtered$description), ]

aldex2_filtered$adjP <- -log10(aldex2_filtered$BH_Corrected_P_Value)

aldex2_filtered_sorted <- aldex2_filtered[order(-aldex2_filtered$BH_Corrected_P_Value), ]

pdf(paste0(Sys.Date(), "_MM_AC_females_picrust2_barplot_down_in_EF.pdf"), width = 15, height = 5)
ggplot(aldex2_filtered_sorted, aes(x = reorder(description, adjP), y = adjP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + 
  labs(x = "description", y = "-log10(adjP)", title = "Pathway vs. Adjusted P-value") +
  theme_minimal()
dev.off()

aldex2_filtered <- subset(merged_df, average_effect <= -0.5)
aldex2_filtered <- aldex2_filtered[!is.na(aldex2_filtered$description), ]

aldex2_filtered$adjP <- -log10(aldex2_filtered$BH_Corrected_P_Value)

aldex2_filtered_sorted <- aldex2_filtered[order(-aldex2_filtered$adjP), ]

pdf(paste0(Sys.Date(), "_MM_AC_females_picrust2_barplot_up_in_EF.pdf"), width = 15, height = 5)
ggplot(aldex2_filtered_sorted, aes(x = reorder(description, adjP), y = adjP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + 
  labs(x = "description", y = "-log10(adjP)", title = "Pathway vs. Adjusted P-value") +
  theme_minimal()
dev.off()

###################################
# 4. Assess GUS activity - estrobolome
###################################

# Access beta-glucuronidase data

AC16.bg <- AC.females.16.pathway.pred["EC:3.2.1.31",]
AC17.bg <- AC.females.17.pathway.pred["EC:3.2.1.31",]
AC22.bg <- AC.females.22.pathway.pred["EC:3.2.1.31",]
AC24.bg <- AC.females.24.pathway.pred["EC:3.2.1.31",]

normalize_by_4m_median <- function(bg_vector, cohort_name) {

  is_4m <- grepl("_4m_", names(bg_vector))
  is_20m <- grepl("_20m_", names(bg_vector))
  
  median_4m <- median(as.numeric(bg_vector[is_4m]), na.rm = TRUE)
  
  norm_values <- as.numeric(bg_vector) / median_4m
  
  data.frame(
    Sample = names(bg_vector),
    Normalized = norm_values,
    Age = ifelse(is_4m, "4m", "20m"),
    Cohort = cohort_name,
    stringsAsFactors = FALSE
  )
}

# Apply function to each cohort
df_16 <- normalize_by_4m_median(AC16.bg, "AC16")
df_17 <- normalize_by_4m_median(AC17.bg, "AC17")
df_22 <- normalize_by_4m_median(AC22.bg, "AC22")
df_24 <- normalize_by_4m_median(AC24.bg, "AC24")

# Combine into one data frame
combined_df <- rbind(df_16, df_17, df_22, df_24)
combined_df$Age <- factor(combined_df$Age, levels = c("4m", "20m"))
combined_df$Cohort <- factor(combined_df$Cohort, levels = c("AC16", "AC17", "AC22", "AC24"))

pch_map <- c("AC16" = 16, "AC17" = 17, "AC22" = 15, "AC24" = 18)
point_shapes <- pch_map[as.character(combined_df$Cohort)]

wilcox.test(combined_df$Normalized[combined_df$Age == "4m"], combined_df$Normalized[combined_df$Age == "20m"])

pdf(paste(Sys.Date(), "MM_AC_cohort_GUS_activity_quantification.pdf", sep = "_"), width = 6, height = 7)
boxplot(Normalized ~ Age, data = combined_df, outline = FALSE,
        main = "GUS activity",
        col = c("deeppink", "deeppink4"),
        xlab = "Age",
        ylab = "GUS activity",
        ylim = c(0.5, 2))
beeswarm(Normalized ~ Age, data = combined_df, pch = 16, add = TRUE)
text(1.5, 2, 
     paste("p-value ~0.3759"))
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_AC_females_differential_abundance_analysis.txt", sep =""))
sessionInfo()
sink()
