setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Justin/Codes/0_Done/2_16S_amplicon_data_analysis/4_Alpha_and_beta_diversity_analysis/0_Generate_plots_R/")

library('qiime2R')
library(tidyverse)

###################################
# Menopause-microbiome project
# Generate beta diversity plots
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

my.braycurtis.AC16.females  <- read_qza("../MM_AC16_females_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")
my.braycurtis.AC17.females  <- read_qza("../MM_AC17_females_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")
my.braycurtis.AC22.females  <- read_qza("../MM_AC22_females_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")
my.braycurtis.AC24.females  <- read_qza("../MM_AC24_females_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")

my.jaccard.AC16.females  <- read_qza("../MM_AC16_females_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")
my.jaccard.AC17.females  <- read_qza("../MM_AC17_females_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")
my.jaccard.AC22.females  <- read_qza("../MM_AC22_females_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")
my.jaccard.AC24.females  <- read_qza("../MM_AC24_females_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")

# 2. Generate plots

my.braycurtis.AC16.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC16) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC16_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.braycurtis.AC17.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC17) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC17_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.braycurtis.AC22.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC22) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC22_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.braycurtis.AC24.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC24) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC24_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 


my.jaccard.AC16.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC16) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC16_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.AC17.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC17) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC17_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.AC22.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC22) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC22_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.AC24.females$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC24) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_females_AC24_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

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

# 2. Generate plots

my.braycurtis.AC16.males$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC16) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_males_AC16_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.braycurtis.AC17.males$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC17) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_males_AC17_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.braycurtis.AC22.males$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC22) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_males_AC22_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.AC16.males$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC16) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_males_AC16_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.AC17.males$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC17) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_males_AC17_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.AC22.males$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.AC22) %>%
        ggplot(aes(x=PC1, y=PC2, color=`age_cat`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="age_cat", values=c("deeppink1", "deeppink4"))
ggsave("MM_AC_males_AC22_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

##############################################
# VCD cohort
##############################################

# 1. Import data files 

my.metadata.VCD <- read_q2metadata("../../1_Manifest_and_metadata/MM_VCD_cohort_sample-metadata.tsv")
my.metadata.VCD$treatment <- factor(my.metadata.VCD$treatment, levels = c("CTL", "VCD"))

my.braycurtis.VCD <- read_qza("../MM_VCD_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")
my.jaccard.VCD    <- read_qza("../MM_VCD_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")

# 2. Generate plots

my.braycurtis.VCD$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.VCD) %>%
        ggplot(aes(x=PC1, y=PC2, color=`treatment`)) +
        geom_point(alpha=0.9, size = 3) + 
        theme_q2r() +
        scale_color_manual(name="treatment", values=c("#7AD151FF", "#FDE725FF"))
ggsave("MM_VCD_PostM_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf") 

my.jaccard.VCD$data$Vectors %>%
        select(SampleID, PC1, PC3) %>%
        left_join(my.metadata.VCD) %>%
        ggplot(aes(x=PC1, y=PC3, color=`treatment`)) +
        geom_point(alpha=0.9, size = 3) +
        theme_q2r() +
        scale_color_manual(name="treatment", values=c("#7AD151FF", "#FDE725FF"))
ggsave("MM_VCD_PostM_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

##############################################
# FMT cohort
##############################################

# 1. Import data files 

my.metadata.FMT <- read_q2metadata("../../1_Manifest_and_metadata/MM_FMT_AC_PostFMT_sample-metadata.tsv")
my.metadata.FMT$fmt <- factor(my.metadata.FMT$fmt, levels = c("FMT_Y", "FMT_O"))

my.braycurtis.FMT <- read_qza("../MM_FMT_AC_PostFMT_diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza")
my.jaccard.FMT    <- read_qza("../MM_FMT_AC_PostFMT_diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")

# 2. Generate plots

my.braycurtis.FMT$data$Vectors %>%
        select(SampleID, PC1, PC2) %>%
        left_join(my.metadata.FMT) %>%
        ggplot(aes(x=PC1, y=PC2, color=`fmt`)) +
        geom_point(alpha=0.9, size = 3) +
        theme_q2r() +
        scale_color_manual(name="fmt", values=c("purple", "orange"))
ggsave("MM_FMT_PostM_Bray_Curtis_index_PCoA.pdf", height=4, width=5, device="pdf")

my.jaccard.FMT$data$Vectors %>%
        select(SampleID, PC1, PC3) %>%
        left_join(my.metadata.FMT) %>%
        ggplot(aes(x=PC1, y=PC3, color=`fmt`)) +
        geom_point(alpha=0.9, size = 3) +
        theme_q2r() +
        scale_color_manual(name="fmt", values=c("purple", "orange"))
ggsave("MM_FMT_PostM_Jaccard_index_PCoA.pdf", height=4, width=5, device="pdf") 

###########################################################
sink(file = paste(Sys.Date(),"MM_beta_diversity_plots.txt", sep =""))
sessionInfo()
sink()
