setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Revision_nat_aging/Data/VCD_combined/4_Alpha_and_beta_diversity_analysis/")

library('qiime2R')
library(tidyverse)
library(mixOmics)
library(PLSDAbatch)

##############################################
# Menopause-microbiome project - revision
# VCD cohort
# https://evayiwenwang.github.io/Managing_batch_effects/index.html#data-processing
##############################################

##############################################
# 1. Import data files
##############################################

# Import files

metadata.VCD <- read_q2metadata("../1_Manifest_and_metadata/MM_VCD_cohort_combined_sample-metadata.tsv")
metadata.VCD$treatment <- factor(metadata.VCD$treatment, levels = c("CTL", "VCD"))

###############################
# 2. Plot Bray-Curtis & Jaccard (beta diversity) plots
###############################

# Create sex_group column
metadata.VCD <- metadata.VCD %>%
  mutate(sex_group = paste(sex, treatment, sep = "_"))

metadata.VCD$sex_group       <- factor(metadata.VCD$sex_group, levels = c("F_CTL", "F_VCD", "M_CTL", "M_VCD"))

my.braycurtis  <- read_qza("./MM_VCD_combined_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")
my.jaccard     <- read_qza("./MM_VCD_combined_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")

# 2. Generate plots

my.braycurtis$data$Vectors %>%
  select(SampleID, PC1, PC3) %>%
  left_join(metadata.VCD) %>%
  ggplot(aes(x=PC1, y=PC3, color=`sex_group`)) +
  geom_point(alpha=0.9, size = 3) + 
  theme_q2r() +
  scale_color_manual(name="sex_group", values=c("#7AD151FF", "#FDE725FF", "#3B9AD8", "#F28E2B"))
ggsave("MM_VCD_combined_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard$data$Vectors %>%
  select(SampleID, PC1, PC3) %>%
  left_join(metadata.VCD) %>%
  ggplot(aes(x=PC1, y=PC3, color=`sex_group`)) +
  geom_point(alpha=0.9, size = 3) + 
  theme_q2r() +
  scale_color_manual(name="sex_group", values=c("#7AD151FF", "#FDE725FF", "#3B9AD8", "#F28E2B"))
ggsave("MM_VCD_combined_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

###############################
# 3. Plot alpha diversity plots
###############################

my.observed.features <- read_qza("~/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Revision_nat_aging/Data/VCD_combined/4_Alpha_and_beta_diversity_analysis_1/MM_VCD_combined_diversity-core-metrics-phylogenetic/observed_features_vector.qza")
my.shannon           <- read_qza("~/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Revision_nat_aging/Data/VCD_combined/4_Alpha_and_beta_diversity_analysis_1/MM_VCD_combined_diversity-core-metrics-phylogenetic/shannon_vector.qza")

my.observed.features <- my.observed.features$data %>% rownames_to_column("SampleID") 
my.shannon           <- my.shannon$data %>% rownames_to_column("SampleID") 

# Normalize data by median

# Calculate median

observed.features.median <- median(my.observed.features$observed_feature)
shannon.median           <- median(my.shannon$shannon_entropy)

my.observed.features$observed_features_normalized <- my.observed.features$observed_features/observed.features.median
my.shannon$shannon_entropy_normalized             <- my.shannon$shannon_entropy/shannon.median

my.alpha.diversity.combined <- cbind(my.observed.features, my.shannon)
rownames(my.alpha.diversity.combined) <- my.alpha.diversity.combined$SampleID

my.alpha.diversity.combined <- my.alpha.diversity.combined[,c(2,3,5,6)]
my.alpha.diversity.combined$Group <- my.metadata$group
my.alpha.diversity.combined$Sex <- my.metadata$sex

my.alpha.diversity.combined.male <- my.alpha.diversity.combined[(my.alpha.diversity.combined$Sex == "M"),]

# Plot data 

observed.features <- list("CTL" = my.alpha.diversity.combined.male$observed_features_normalized[grep("CTL", my.alpha.diversity.combined.male$Group)],
                          "VCD" = my.alpha.diversity.combined.male$observed_features_normalized[grep("VCD", my.alpha.diversity.combined.male$Group)])

shannon.entropy   <- list("CTL" = my.alpha.diversity.combined.male$shannon_entropy_normalized[grep("CTL", my.alpha.diversity.combined.male$Group)],
                          "VCD" = my.alpha.diversity.combined.male$shannon_entropy_normalized[grep("VCD", my.alpha.diversity.combined.male$Group)])


wilcox.test(my.alpha.diversity.combined.male$observed_features_normalized[grep("CTL", my.alpha.diversity.combined.male$Group)], 
            my.alpha.diversity.combined.male$observed_features_normalized[grep("VCD", my.alpha.diversity.combined.male$Group)],)       # p-value = 0.4286
wilcox.test(my.alpha.diversity.combined.male$shannon_entropy_normalized[grep("CTL", my.alpha.diversity.combined.male$Group)], 
            my.alpha.diversity.combined.male$shannon_entropy_normalized[grep("VCD", my.alpha.diversity.combined.male$Group)],)         # p-value = 0.4286

pdf(paste0(Sys.Date(), "_MM_VCD_Cohort_combined_normalized_alpha_diversity_observed_features.pdf"), width = 4, height = 7)
boxplot(observed.features, 
        outline = TRUE, 
        col = c("blue", "red"),
        main = "Observed features (normalized)", ylab = "Observed features (normalized)", xlab = "Age", ylim = c(0.5, 1.5))
beeswarm(observed.features, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.4286", cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_MM_VCD_Cohort_combined_normalized_alpha_diversity_shannon_entropy.pdf"), width = 4, height = 7)
boxplot(shannon.entropy, 
        outline = TRUE, 
        col = c("blue", "red"),
        main = "Shannon entropy (normalized)", ylab = "Shannon entropy (normalized)", xlab = "Age",
        ylim = c(0.5, 1.5))
beeswarm(shannon.entropy, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5,1.5,"p-value ~0.4286", cex = 1)
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_VCD_diversity_plots.txt", sep =""))
sessionInfo()
sink()
