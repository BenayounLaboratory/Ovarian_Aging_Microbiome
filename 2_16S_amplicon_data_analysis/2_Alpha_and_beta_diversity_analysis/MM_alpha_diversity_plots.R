setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Github/16S_amplicon/4_Alpha_and_beta_diversity_analysis/0_Generate_plots_R/")

library('qiime2R')
library(tidyverse)
library('beeswarm')
# library(ggrepel) # for offset labels
# library(ggtree) # for visualizing phylogenetic trees
# library(ape) # for manipulating phylogenetic trees

###################################
# Menopause-microbiome project
# Generate alpha diversity plots
###################################

##############################################
# AC cohort - females
##############################################

# 1. Import data files 

my.metadata.AC16 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC16_AC_females_sample-metadata.txt")
my.metadata.AC17 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC17_AC_females_sample-metadata.txt")
my.metadata.AC22 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC22_AC_females_sample-metadata.txt")
my.metadata.AC24 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC24_AC_females_sample-metadata.txt")

my.metadata.AC16$age_cat <- factor(my.metadata.AC16$age_cat, levels = c("4m", "20m"))
my.metadata.AC17$age_cat <- factor(my.metadata.AC17$age_cat, levels = c("4m", "20m"))
my.metadata.AC22$age_cat <- factor(my.metadata.AC22$age_cat, levels = c("4m", "20m"))
my.metadata.AC24$age_cat <- factor(my.metadata.AC24$age_cat, levels = c("4m", "20m"))

my.metadata.combined <- rbind(my.metadata.AC16, my.metadata.AC17, my.metadata.AC22, my.metadata.AC24)

my.SVs.AC16 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC16_females_table.qza")$data
my.SVs.AC17 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC17_females_table.qza")$data
my.SVs.AC22 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC22_females_table.qza")$data
my.SVs.AC24 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC24_females_table.qza")$data

my.observed.features.AC16 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_females_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC16           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_females_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC16          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_females_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC16           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_females_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC17 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_females_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC17           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_females_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC17          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_females_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC17           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_females_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC22 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_females_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC22           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_females_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC22          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_females_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC22           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_females_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC24 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC24_females_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC24           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC24_females_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC24          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC24_females_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC24           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC24_females_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC16 <- my.observed.features.AC16$data %>% rownames_to_column("SampleID") 
my.shannon.AC16           <- my.shannon.AC16$data %>% rownames_to_column("SampleID") 
my.evenness.AC16          <- my.evenness.AC16$data %>% rownames_to_column("SampleID") 
my.faithpd.AC16           <- my.faithpd.AC16$data

my.observed.features.AC17 <- my.observed.features.AC17$data %>% rownames_to_column("SampleID") 
my.shannon.AC17           <- my.shannon.AC17$data %>% rownames_to_column("SampleID") 
my.evenness.AC17          <- my.evenness.AC17$data %>% rownames_to_column("SampleID") 
my.faithpd.AC17           <- my.faithpd.AC17$data

my.observed.features.AC22 <- my.observed.features.AC22$data %>% rownames_to_column("SampleID") 
my.shannon.AC22           <- my.shannon.AC22$data %>% rownames_to_column("SampleID") 
my.evenness.AC22          <- my.evenness.AC22$data %>% rownames_to_column("SampleID") 
my.faithpd.AC22           <- my.faithpd.AC22$data

my.observed.features.AC24 <- my.observed.features.AC24$data %>% rownames_to_column("SampleID") 
my.shannon.AC24           <- my.shannon.AC24$data %>% rownames_to_column("SampleID") 
my.evenness.AC24          <- my.evenness.AC24$data %>% rownames_to_column("SampleID") 
my.faithpd.AC24           <- my.faithpd.AC24$data

colnames(my.faithpd.AC16) <- c("SampleID", "faithpd")
colnames(my.faithpd.AC17) <- c("SampleID", "faithpd")
colnames(my.faithpd.AC22) <- c("SampleID", "faithpd")
colnames(my.faithpd.AC24) <- c("SampleID", "faithpd")

# 3. Normalize data by median

# Calculate median

# Cohort 16

observed.features.median.AC16 <- median(my.observed.features.AC16$observed_feature)
shannon.median.AC16           <- median(my.shannon.AC16$shannon_entropy)
evenness.median.AC16          <- median(my.evenness.AC16$pielou_evenness)
faithpd.median.AC16           <- median(my.faithpd.AC16$faithpd)

# Cohort 17

observed.features.median.AC17 <- median(my.observed.features.AC17$observed_feature)
shannon.median.AC17           <- median(my.shannon.AC17$shannon_entropy)
evenness.median.AC17          <- median(my.evenness.AC17$pielou_evenness)
faithpd.median.AC17           <- median(my.faithpd.AC17$faithpd)

# Cohort 22

observed.features.median.AC22 <- median(my.observed.features.AC22$observed_feature)
shannon.median.AC22           <- median(my.shannon.AC22$shannon_entropy)
evenness.median.AC22          <- median(my.evenness.AC22$pielou_evenness)
faithpd.median.AC22           <- median(my.faithpd.AC22$faithpd)

# Cohort 24

observed.features.median.AC24 <- median(my.observed.features.AC24$observed_feature)
shannon.median.AC24           <- median(my.shannon.AC24$shannon_entropy)
evenness.median.AC24          <- median(my.evenness.AC24$pielou_evenness)
faithpd.median.AC24           <- median(my.faithpd.AC24$faithpd)

# Normalize by median

# Cohort 16

my.observed.features.AC16$observed_features_normalized <- my.observed.features.AC16$observed_features/observed.features.median.AC16
my.shannon.AC16$shannon_entropy_normalized             <- my.shannon.AC16$shannon_entropy/shannon.median.AC16
my.evenness.AC16$pielou_evenness_normalized            <- my.evenness.AC16$pielou_evenness/evenness.median.AC16
my.faithpd.AC16$faithpd_normalized                     <- my.faithpd.AC16$faithpd/faithpd.median.AC16

# Cohort 17

my.observed.features.AC17$observed_features_normalized <- my.observed.features.AC17$observed_features/observed.features.median.AC17
my.shannon.AC17$shannon_entropy_normalized             <- my.shannon.AC17$shannon_entropy/shannon.median.AC17
my.evenness.AC17$pielou_evenness_normalized            <- my.evenness.AC17$pielou_evenness/evenness.median.AC17
my.faithpd.AC17$faithpd_normalized                     <- my.faithpd.AC17$faithpd/faithpd.median.AC17

# Cohort 22

my.observed.features.AC22$observed_features_normalized <- my.observed.features.AC22$observed_features/observed.features.median.AC22
my.shannon.AC22$shannon_entropy_normalized             <- my.shannon.AC22$shannon_entropy/shannon.median.AC22
my.evenness.AC22$pielou_evenness_normalized            <- my.evenness.AC22$pielou_evenness/evenness.median.AC22
my.faithpd.AC22$faithpd_normalized                     <- my.faithpd.AC22$faithpd/faithpd.median.AC22

# Cohort 24

my.observed.features.AC24$observed_features_normalized <- my.observed.features.AC24$observed_features/observed.features.median.AC24
my.shannon.AC24$shannon_entropy_normalized             <- my.shannon.AC24$shannon_entropy/shannon.median.AC24
my.evenness.AC24$pielou_evenness_normalized            <- my.evenness.AC24$pielou_evenness/evenness.median.AC24
my.faithpd.AC24$faithpd_normalized                     <- my.faithpd.AC24$faithpd/faithpd.median.AC24

# Combine data

my.alpha.diversity.AC16 <- cbind(my.observed.features.AC16, my.shannon.AC16[,2:3], my.evenness.AC16[,2:3], my.faithpd.AC16[,2:3])
my.alpha.diversity.AC17 <- cbind(my.observed.features.AC17, my.shannon.AC17[,2:3], my.evenness.AC17[,2:3], my.faithpd.AC17[,2:3])
my.alpha.diversity.AC22 <- cbind(my.observed.features.AC22, my.shannon.AC22[,2:3], my.evenness.AC22[,2:3], my.faithpd.AC22[,2:3])
my.alpha.diversity.AC24 <- cbind(my.observed.features.AC24, my.shannon.AC24[,2:3], my.evenness.AC24[,2:3], my.faithpd.AC24[,2:3])

my.alpha.diversity.combined <- rbind(my.alpha.diversity.AC16, my.alpha.diversity.AC17, my.alpha.diversity.AC22, my.alpha.diversity.AC24)
rownames(my.alpha.diversity.combined) <- my.alpha.diversity.combined$SampleID

# 4. Plot data 

observed.features <- list("YF" = my.alpha.diversity.combined$observed_features_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "EF" = my.alpha.diversity.combined$observed_features_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

shannon.entropy   <- list("YF" = my.alpha.diversity.combined$shannon_entropy_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "EF" = my.alpha.diversity.combined$shannon_entropy_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

evenness          <- list("YF" = my.alpha.diversity.combined$pielou_evenness_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "EF" = my.alpha.diversity.combined$pielou_evenness_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

faithpd           <- list("YF" = my.alpha.diversity.combined$faithpd_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "EF" = my.alpha.diversity.combined$faithpd_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

wilcox.test(my.alpha.diversity.combined$observed_features_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$observed_features_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 0.0009126
wilcox.test(my.alpha.diversity.combined$shannon_entropy_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$shannon_entropy_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 9.652e-06
wilcox.test(my.alpha.diversity.combined$pielou_evenness_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$pielou_evenness_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 2.604e-05
wilcox.test(my.alpha.diversity.combined$faithpd_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$faithpd_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 0.01387

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_females_combined_normalized_alpha_diversity_observed_features.pdf"), width = 4, height = 7)
boxplot(observed.features, 
        outline = TRUE, 
        col = c("deeppink1", "deeppink4"),
        main = "Observed features (normalized)", ylab = "Observed features (normalized)", xlab = "Age", ylim = c(0.5, 1.5))
beeswarm(observed.features, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.0009126", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_females_combined_normalized_alpha_diversity_shannon_entropy.pdf"), width = 4, height = 7)
boxplot(shannon.entropy, 
        outline = TRUE, 
        col = c("deeppink1", "deeppink4"),
        main = "Shannon entropy (normalized)", ylab = "Shannon entropy (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(shannon.entropy, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~9.652e-06", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_females_combined_normalized_alpha_diversity_evenness_pielou.pdf"), width = 4, height = 7)
boxplot(evenness, 
        outline = TRUE, 
        col = c("deeppink1", "deeppink4"),
        main = "Evenness pielou (normalized)", ylab = "Evenness pielou (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(evenness, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~2.604e-05", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_females_combined_normalized_alpha_diversity_faithpd.pdf"), width = 4, height = 7)
boxplot(faithpd, 
        outline = TRUE, 
        col = c("deeppink1", "deeppink4"),
        main = "Faith PD (normalized)", ylab = "Faith PD (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(faithpd, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 0.01387", cex = 1)
dev.off()

##############################################
# AC cohort - males
##############################################

# 1. Import data files 

my.metadata.AC16 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC16_AC_males_sample-metadata.txt")
my.metadata.AC17 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC17_AC_males_sample-metadata.txt")
my.metadata.AC22 <- read_q2metadata("../../1_Manifest_and_metadata/MM_AC22_AC_males_sample-metadata.txt")

my.metadata.AC16$age_cat <- factor(my.metadata.AC16$age_cat, levels = c("4m", "20m"))
my.metadata.AC17$age_cat <- factor(my.metadata.AC17$age_cat, levels = c("4m", "20m"))
my.metadata.AC22$age_cat <- factor(my.metadata.AC22$age_cat, levels = c("4m", "20m"))

my.metadata.combined <- rbind(my.metadata.AC16, my.metadata.AC17, my.metadata.AC22)

my.SVs.AC16 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC16_males_table.qza")$data
my.SVs.AC17 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC17_males_table.qza")$data
my.SVs.AC22 <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_AC22_males_table.qza")$data

my.observed.features.AC16 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_males_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC16           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_males_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC16          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_males_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC16           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC16_males_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC17 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_males_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC17           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_males_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC17          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_males_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC17           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC17_males_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC22 <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_males_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.AC22           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_males_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.AC22          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_males_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.AC22           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_AC22_males_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.AC16 <- my.observed.features.AC16$data %>% rownames_to_column("SampleID") 
my.shannon.AC16           <- my.shannon.AC16$data %>% rownames_to_column("SampleID") 
my.evenness.AC16          <- my.evenness.AC16$data %>% rownames_to_column("SampleID") 
my.faithpd.AC16           <- my.faithpd.AC16$data

my.observed.features.AC17 <- my.observed.features.AC17$data %>% rownames_to_column("SampleID") 
my.shannon.AC17           <- my.shannon.AC17$data %>% rownames_to_column("SampleID") 
my.evenness.AC17          <- my.evenness.AC17$data %>% rownames_to_column("SampleID") 
my.faithpd.AC17           <- my.faithpd.AC17$data

my.observed.features.AC22 <- my.observed.features.AC22$data %>% rownames_to_column("SampleID") 
my.shannon.AC22           <- my.shannon.AC22$data %>% rownames_to_column("SampleID") 
my.evenness.AC22          <- my.evenness.AC22$data %>% rownames_to_column("SampleID") 
my.faithpd.AC22           <- my.faithpd.AC22$data

colnames(my.faithpd.AC16) <- c("SampleID", "faithpd")
colnames(my.faithpd.AC17) <- c("SampleID", "faithpd")
colnames(my.faithpd.AC22) <- c("SampleID", "faithpd")

# 3. Normalize data by median

# Calculate median

# Cohort 16

observed.features.median.AC16 <- median(my.observed.features.AC16$observed_feature)
shannon.median.AC16           <- median(my.shannon.AC16$shannon_entropy)
evenness.median.AC16          <- median(my.evenness.AC16$pielou_evenness)
faithpd.median.AC16           <- median(my.faithpd.AC16$faithpd)

# Cohort 17

observed.features.median.AC17 <- median(my.observed.features.AC17$observed_feature)
shannon.median.AC17           <- median(my.shannon.AC17$shannon_entropy)
evenness.median.AC17          <- median(my.evenness.AC17$pielou_evenness)
faithpd.median.AC17           <- median(my.faithpd.AC17$faithpd)

# Cohort 22

observed.features.median.AC22 <- median(my.observed.features.AC22$observed_feature)
shannon.median.AC22           <- median(my.shannon.AC22$shannon_entropy)
evenness.median.AC22          <- median(my.evenness.AC22$pielou_evenness)
faithpd.median.AC22           <- median(my.faithpd.AC22$faithpd)

# Normalize by median

# Cohort 16

my.observed.features.AC16$observed_features_normalized <- my.observed.features.AC16$observed_features/observed.features.median.AC16
my.shannon.AC16$shannon_entropy_normalized             <- my.shannon.AC16$shannon_entropy/shannon.median.AC16
my.evenness.AC16$pielou_evenness_normalized            <- my.evenness.AC16$pielou_evenness/evenness.median.AC16
my.faithpd.AC16$faithpd_normalized                     <- my.faithpd.AC16$faithpd/faithpd.median.AC16

# Cohort 17

my.observed.features.AC17$observed_features_normalized <- my.observed.features.AC17$observed_features/observed.features.median.AC17
my.shannon.AC17$shannon_entropy_normalized             <- my.shannon.AC17$shannon_entropy/shannon.median.AC17
my.evenness.AC17$pielou_evenness_normalized            <- my.evenness.AC17$pielou_evenness/evenness.median.AC17
my.faithpd.AC17$faithpd_normalized                     <- my.faithpd.AC17$faithpd/faithpd.median.AC17

# Cohort 22

my.observed.features.AC22$observed_features_normalized <- my.observed.features.AC22$observed_features/observed.features.median.AC22
my.shannon.AC22$shannon_entropy_normalized             <- my.shannon.AC22$shannon_entropy/shannon.median.AC22
my.evenness.AC22$pielou_evenness_normalized            <- my.evenness.AC22$pielou_evenness/evenness.median.AC22
my.faithpd.AC22$faithpd_normalized                     <- my.faithpd.AC22$faithpd/faithpd.median.AC22

# Combine data

my.alpha.diversity.AC16 <- cbind(my.observed.features.AC16, my.shannon.AC16[,2:3], my.evenness.AC16[,2:3], my.faithpd.AC16[,2:3])
my.alpha.diversity.AC17 <- cbind(my.observed.features.AC17, my.shannon.AC17[,2:3], my.evenness.AC17[,2:3], my.faithpd.AC17[,2:3])
my.alpha.diversity.AC22 <- cbind(my.observed.features.AC22, my.shannon.AC22[,2:3], my.evenness.AC22[,2:3], my.faithpd.AC22[,2:3])

my.alpha.diversity.combined <- rbind(my.alpha.diversity.AC16, my.alpha.diversity.AC17, my.alpha.diversity.AC22)
rownames(my.alpha.diversity.combined) <- my.alpha.diversity.combined$SampleID

# 4. Plot data 

observed.features <- list("YM" = my.alpha.diversity.combined$observed_features_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "OM" = my.alpha.diversity.combined$observed_features_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

shannon.entropy   <- list("YM" = my.alpha.diversity.combined$shannon_entropy_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "OM" = my.alpha.diversity.combined$shannon_entropy_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

evenness          <- list("YM" = my.alpha.diversity.combined$pielou_evenness_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "OM" = my.alpha.diversity.combined$pielou_evenness_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

faithpd           <- list("YM" = my.alpha.diversity.combined$faithpd_normalized[grep("4m", my.alpha.diversity.combined$SampleID)],
                          "OM" = my.alpha.diversity.combined$faithpd_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])

wilcox.test(my.alpha.diversity.combined$observed_features_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$observed_features_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 0.9633
wilcox.test(my.alpha.diversity.combined$shannon_entropy_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$shannon_entropy_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 0.2456
wilcox.test(my.alpha.diversity.combined$pielou_evenness_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$pielou_evenness_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 0.1499
wilcox.test(my.alpha.diversity.combined$faithpd_normalized[grep("4m", my.alpha.diversity.combined$SampleID)], my.alpha.diversity.combined$faithpd_normalized[grep("20m", my.alpha.diversity.combined$SampleID)])     # p-value = 0.8743

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_males_combined_normalized_alpha_diversity_observed_features.pdf"), width = 4, height = 7)
boxplot(observed.features, 
        outline = TRUE, 
        col = c("deepskyblue", "deepskyblue4"),
        main = "Observed features (normalized)", ylab = "Observed features (normalized)", xlab = "Age", ylim = c(0.5, 1.5))
beeswarm(observed.features, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.9633", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_males_combined_normalized_alpha_diversity_shannon_entropy.pdf"), width = 4, height = 7)
boxplot(shannon.entropy, 
        outline = TRUE, 
        col = c("deepskyblue", "deepskyblue4"),
        main = "Shannon entropy (normalized)", ylab = "Shannon entropy (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(shannon.entropy, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.2456", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_males_combined_normalized_alpha_diversity_evenness_pielou.pdf"), width = 4, height = 7)
boxplot(evenness, 
        outline = TRUE, 
        col = c("deepskyblue", "deepskyblue4"),
        main = "Evenness pielou (normalized)", ylab = "Evenness pielou (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(evenness, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.1499", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_AC_Cohorts_males_combined_normalized_alpha_diversity_faithpd.pdf"), width = 4, height = 7)
boxplot(faithpd, 
        outline = TRUE, 
        col = c("deepskyblue", "deepskyblue4"),
        main = "Faith PD (normalized)", ylab = "Faith PD (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(faithpd, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 0.8743", cex = 1)
dev.off()

##############################################
# VCD cohort
##############################################

# 1. Import data files 

my.metadata.VCD <- read_q2metadata("../../1_Manifest_and_metadata/MM_VCD_cohort_sample-metadata.tsv")

my.metadata.VCD$treatment <- factor(my.metadata.VCD$treatment, levels = c("CTL", "VCD"))

my.SVs.VCD <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_VCD_table.qza")$data

my.observed.features.VCD <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_VCD_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.VCD           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_VCD_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.VCD          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_VCD_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.VCD           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_VCD_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.VCD <- my.observed.features.VCD$data %>% rownames_to_column("SampleID") 
my.shannon.VCD           <- my.shannon.VCD$data %>% rownames_to_column("SampleID") 
my.evenness.VCD          <- my.evenness.VCD$data %>% rownames_to_column("SampleID") 
my.faithpd.VCD           <- my.faithpd.VCD$data

colnames(my.faithpd.VCD) <- c("SampleID", "faithpd")

# 3. Normalize data by median

# Calculate median

observed.features.median.VCD <- median(my.observed.features.VCD$observed_feature)
shannon.median.VCD           <- median(my.shannon.VCD$shannon_entropy)
evenness.median.VCD          <- median(my.evenness.VCD$pielou_evenness)
faithpd.median.VCD           <- median(my.faithpd.VCD$faithpd)

# Normalize by median

my.observed.features.VCD$observed_features_normalized <- my.observed.features.VCD$observed_features/observed.features.median.VCD
my.shannon.VCD$shannon_entropy_normalized             <- my.shannon.VCD$shannon_entropy/shannon.median.VCD
my.evenness.VCD$pielou_evenness_normalized            <- my.evenness.VCD$pielou_evenness/evenness.median.VCD
my.faithpd.VCD$faithpd_normalized                     <- my.faithpd.VCD$faithpd/faithpd.median.VCD

# Combine data

my.alpha.diversity.VCD <- cbind(my.metadata.VCD$treatment, my.observed.features.VCD, my.shannon.VCD[,2:3], my.evenness.VCD[,2:3], my.faithpd.VCD[,2:3])

# 4. Plot data 

observed.features <- list("CTL" = my.alpha.diversity.VCD$observed_features_normalized[grep("CTL", my.alpha.diversity.VCD[,1])],
                          "VCD" = my.alpha.diversity.VCD$observed_features_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])

shannon.entropy   <- list("CTL" = my.alpha.diversity.VCD$shannon_entropy_normalized[grep("CTL", my.alpha.diversity.VCD[,1])],
                          "VCD" = my.alpha.diversity.VCD$shannon_entropy_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])

evenness          <- list("CTL" = my.alpha.diversity.VCD$pielou_evenness_normalized[grep("CTL", my.alpha.diversity.VCD[,1])],
                          "VCD" = my.alpha.diversity.VCD$pielou_evenness_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])

faithpd           <- list("CTL" = my.alpha.diversity.VCD$faithpd_normalized[grep("CTL", my.alpha.diversity.VCD[,1])],
                          "VCD" = my.alpha.diversity.VCD$faithpd_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])

wilcox.test(my.alpha.diversity.VCD$observed_features_normalized[grep("CTL", my.alpha.diversity.VCD[,1])], my.alpha.diversity.VCD$observed_features_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])     # p-value = 0.6905
wilcox.test(my.alpha.diversity.VCD$shannon_entropy_normalized[grep("CTL", my.alpha.diversity.VCD[,1])], my.alpha.diversity.VCD$shannon_entropy_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])     # p-value = 1
wilcox.test(my.alpha.diversity.VCD$pielou_evenness_normalized[grep("CTL", my.alpha.diversity.VCD[,1])], my.alpha.diversity.VCD$pielou_evenness_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])     # p-value = 1
wilcox.test(my.alpha.diversity.VCD$faithpd_normalized[grep("CTL", my.alpha.diversity.VCD[,1])], my.alpha.diversity.VCD$faithpd_normalized[grep("VCD", my.alpha.diversity.VCD[,1])])     # p-value = 1

pdf(paste0(Sys.Date(), "_MM_VCD_combined_normalized_alpha_diversity_observed_features.pdf"), width = 4, height = 7)
boxplot(observed.features, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Observed features (normalized)", ylab = "Observed features (normalized)", xlab = "Age", ylim = c(0.5, 1.5))
beeswarm(observed.features, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.6905", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_normalized_alpha_diversity_shannon_entropy.pdf"), width = 4, height = 7)
boxplot(shannon.entropy, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Shannon entropy (normalized)", ylab = "Shannon entropy (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(shannon.entropy, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 1", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_normalized_alpha_diversity_evenness_pielou.pdf"), width = 4, height = 7)
boxplot(evenness, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Evenness pielou (normalized)", ylab = "Evenness pielou (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(evenness, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 1", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_normalized_alpha_diversity_faithpd.pdf"), width = 4, height = 7)
boxplot(faithpd, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Faith PD (normalized)", ylab = "Faith PD (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(faithpd, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 1", cex = 1)
dev.off()

##############################################
# FMT cohort
##############################################

# 1. Import data files 

my.metadata.FMT <- read_q2metadata("../../1_Manifest_and_metadata/MM_FMT_AC_PostFMT_sample-metadata.tsv")

my.metadata.FMT$fmt <- factor(my.metadata.FMT$fmt, levels = c("FMT_Y", "FMT_O"))

my.SVs.FMT <- read_qza("../../3_Feature_table_and_rep_seq_table/MM_FMT_AC_PostFMT_table.qza")$data

my.observed.features.FMT <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_FMT_AC_PostFMT_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon.FMT           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_FMT_AC_PostFMT_diversity-core-metrics-phylogenetic/shannon_vector.qza")
my.evenness.FMT          <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_FMT_AC_PostFMT_diversity-core-metrics-phylogenetic/evenness_vector.qza")
my.faithpd.FMT           <- read_qza("../../4_Alpha_and_beta_diversity_analysis/MM_FMT_AC_PostFMT_diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

my.observed.features.FMT <- my.observed.features.FMT$data %>% rownames_to_column("SampleID") 
my.shannon.FMT           <- my.shannon.FMT$data %>% rownames_to_column("SampleID") 
my.evenness.FMT          <- my.evenness.FMT$data %>% rownames_to_column("SampleID") 
my.faithpd.FMT           <- my.faithpd.FMT$data

colnames(my.faithpd.FMT) <- c("SampleID", "faithpd")

# 3. Normalize data by median

# Calculate median

observed.features.median.FMT <- median(my.observed.features.FMT$observed_feature)
shannon.median.FMT           <- median(my.shannon.FMT$shannon_entropy)
evenness.median.FMT          <- median(my.evenness.FMT$pielou_evenness)
faithpd.median.FMT           <- median(my.faithpd.FMT$faithpd)

# Normalize by median

my.observed.features.FMT$observed_features_normalized <- my.observed.features.FMT$observed_features/observed.features.median.FMT
my.shannon.FMT$shannon_entropy_normalized             <- my.shannon.FMT$shannon_entropy/shannon.median.FMT
my.evenness.FMT$pielou_evenness_normalized            <- my.evenness.FMT$pielou_evenness/evenness.median.FMT
my.faithpd.FMT$faithpd_normalized                     <- my.faithpd.FMT$faithpd/faithpd.median.FMT

# Combine data

my.alpha.diversity.FMT <- cbind(my.metadata.FMT$fmt, my.observed.features.FMT, my.shannon.FMT[,2:3], my.evenness.FMT[,2:3], my.faithpd.FMT[,2:3])

# 4. Plot data 

observed.features <- list("FMT_YF" = my.alpha.diversity.FMT$observed_features_normalized[grep("FMT_Y", my.alpha.diversity.FMT[,1])],
                          "FMT_EF" = my.alpha.diversity.FMT$observed_features_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])

shannon.entropy   <- list("FMT_YF" = my.alpha.diversity.FMT$shannon_entropy_normalized[grep("FMT_YF", my.alpha.diversity.FMT[,1])],
                          "FMT_EF" = my.alpha.diversity.FMT$shannon_entropy_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])

evenness          <- list("FMT_YF" = my.alpha.diversity.FMT$pielou_evenness_normalized[grep("FMT_YF", my.alpha.diversity.FMT[,1])],
                          "FMT_EF" = my.alpha.diversity.FMT$pielou_evenness_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])

faithpd           <- list("FMT_YF" = my.alpha.diversity.FMT$faithpd_normalized[grep("FMT_YF", my.alpha.diversity.FMT[,1])],
                          "FMT_EF" = my.alpha.diversity.FMT$faithpd_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])

wilcox.test(my.alpha.diversity.FMT$observed_features_normalized[grep("FMT_Y", my.alpha.diversity.FMT[,1])], my.alpha.diversity.FMT$observed_features_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])     # p-value = 0.8981
wilcox.test(my.alpha.diversity.FMT$shannon_entropy_normalized[grep("FMT_Y", my.alpha.diversity.FMT[,1])], my.alpha.diversity.FMT$shannon_entropy_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])     # p-value = 1
wilcox.test(my.alpha.diversity.FMT$pielou_evenness_normalized[grep("FMT_Y", my.alpha.diversity.FMT[,1])], my.alpha.diversity.FMT$pielou_evenness_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])     # p-value = 1
wilcox.test(my.alpha.diversity.FMT$faithpd_normalized[grep("FMT_Y", my.alpha.diversity.FMT[,1])], my.alpha.diversity.FMT$faithpd_normalized[grep("FMT_O", my.alpha.diversity.FMT[,1])])     # p-value = 0.535

pdf(paste0(Sys.Date(), "_MM_VCD_combined_normalized_alpha_diversity_observed_features.pdf"), width = 4, height = 7)
boxplot(observed.features, 
        outline = TRUE, 
        col = c("purple", "orange"),
        main = "Observed features (normalized)", ylab = "Observed features (normalized)", xlab = "Age", ylim = c(0.5, 1.5))
beeswarm(observed.features, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.8981", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_normalized_alpha_diversity_shannon_entropy.pdf"), width = 4, height = 7)
boxplot(shannon.entropy, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Shannon entropy (normalized)", ylab = "Shannon entropy (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(shannon.entropy, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 1", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_normalized_alpha_diversity_evenness_pielou.pdf"), width = 4, height = 7)
boxplot(evenness, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Evenness pielou (normalized)", ylab = "Evenness pielou (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(evenness, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value = 1", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_normalized_alpha_diversity_faithpd.pdf"), width = 4, height = 7)
boxplot(faithpd, 
        outline = TRUE, 
        col = c("green", "yellow"),
        main = "Faith PD (normalized)", ylab = "Faith PD (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(faithpd, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.535", cex = 1)
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_alpha_diversity_combined_plots.txt", sep =""))
sessionInfo()
sink()
