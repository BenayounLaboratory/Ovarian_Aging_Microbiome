library('mediation')
library('qiime2R')
library(tidyr)
library(dplyr)
library(ggplot2)

######################################
# Mediation analysis:
#                     - Mediators: features (from FMT amplicon data)
#                     - Treatment: FMT_YF vs. FMT_EF
#                     - Outcome: FMT ovarian bulk RNAseq
# https://cran.r-project.org/web/packages/mediation/mediation.pdf
######################################

######################################
# 1. Import data
######################################

# Amplicon data

my.FMT.16S.table         <- read_qza("../3_Feature_table_and_rep_seq_table/MM_FMT_AC_PostFMT_table.qza")$data
my.FMT.16S.differentials <- read_qza("../6_Differential_abundance_and_functional_prediction_analysis/Differential_abundance_aldex2/aldex2_FMT/differentials.qza")$data
my.FMT.16S.taxonomy      <- read_qza("../5_Taxonomy/silva_MM_FMT_AC_PostFMT_taxonomy.qza")$data
my.FMT.16S.metadata      <- read_q2metadata("../1_Manifest_and_metadata/MM_FMT_AC_PostFMT_sample-metadata.tsv")

my.FMT.16S.metadata$fmt <- factor(my.FMT.16S.metadata$fmt, levels = c("FMT_Y", "FMT_O"))

# Significant features
my.sig.features <- read.table("./FMT_differential_abundance_analysis_significant_features.txt", header = TRUE)

# Bulk ovarian RNA-seq

load("../../3_Ovarian_RNAseq_analysis/FMT_cohort/DESeq2/2024-03-29_MM_FMT_Ovaries_DESeq2_Analysis_STAR_post_SVA_normalized_counts.RData")

######################################
# 2. Pre-process data
######################################

# Filter significant features and perform CLR transformation

my.FMT.16S.table.sig <- my.FMT.16S.table[my.sig.features$Feature.ID,]
my.FMT.clr <- apply(log2(my.FMT.16S.table.sig+0.5), 2, function(x) x-mean(x))

# Treatment variable
treatment <- rep(c("FMT_YF", "FMT_OF"), each=7)

# Mediator variable
t.my.FMT.16S.table <- t(my.FMT.clr)

# Re-order rows to match
rownames(t.my.FMT.16S.table) <- c("FMT_YF_1", "FMT_YF_2", 
                                  "FMT_OF_1", "FMT_OF_2", "FMT_OF_3", 
                                  "FMT_YF_3", "FMT_YF_4", "FMT_YF_5", "FMT_YF_6", "FMT_YF_7",
                                  "FMT_OF_4", "FMT_OF_5", "FMT_OF_6", "FMT_OF_7")

rownames <- rownames(t.my.FMT.16S.table)

split_names <- strsplit(rownames, "_")
groups <- sapply(split_names, function(x) x[2])
numbers <- as.numeric(sapply(split_names, function(x) x[3]))

df <- data.frame(rownames, groups, numbers, stringsAsFactors = FALSE)

df_yf <- df[df$groups == "YF",]
df_of <- df[df$groups == "OF",]

df_yf_ordered <- df_yf[order(df_yf$numbers),]
df_of_ordered <- df_of[order(df_of$numbers),]

df_combined_ordered <- rbind(df_yf_ordered, df_of_ordered)

ordered_rownames <- df_combined_ordered$rownames

t.my.FMT.16S.table <- t.my.FMT.16S.table[ordered_rownames, ]

######################################
# 3. Perform mediation analysis
######################################

mediator <- t.my.FMT.16S.table

# Outcome variable - Use PCA 1 coordinates

mds.result <- cmdscale(1-cor(tissue.cts, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
outcome <- mds.result[, 1]

# Feature.ID: d4840c45eb37593108d020fb6a63f8e2

myData <- data.frame(mediator = mediator[,"d4840c45eb37593108d020fb6a63f8e2"],
                     treatment = treatment,
                     outcome = outcome)

# Model for the mediator
mediator.model <- lm(mediator ~ treatment, data = myData)

# Model for the outcome
outcome.model <- lm(outcome ~ treatment + mediator, data = myData)

med.out <- mediate(mediator.model, outcome.model, treat="treatment", mediator="mediator")

summary(med.out)
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.016051     0.018695         0.01  <2e-16 ***
# ADE             0.006692     0.000024         0.01    0.05 *  
# Total Effect    0.022413     0.029453         0.02  <2e-16 ***

results <- data.frame(EffectType = c("ACME", "ADE", "Total Effect"),
                      Estimate = c(0.016051, 0.006692, 0.022413),
                      LowerCI = c(0.018695, 0.000024, 0.029453),
                      UpperCI = c(0.01, 0.01, 0.02))

pdf("./MM_FMT_RNAseq_mediation_test_D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Lachnospiraceae.pdf", width = 8, height = 6)
ggplot(results, aes(x = EffectType, y = Estimate, fill = EffectType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(.9)) +
  theme_minimal() +
  labs(title = "Causal Mediation Analysis Results - D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Lachnospiraceae",
       x = "Effect Type",
       y = "Effect Estimate") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank())
dev.off()

# Feature.ID: e4f32cdc3407edb8423447421078fdae

myData <- data.frame(mediator = mediator[,"e4f32cdc3407edb8423447421078fdae"],
                     treatment = treatment,
                     outcome = outcome)

# Model for the mediator
mediator.model <- lm(mediator ~ treatment, data = myData)

# Model for the outcome
outcome.model <- lm(outcome ~ treatment + mediator, data = myData)

med.out <- mediate(mediator.model, outcome.model, treat="treatment", mediator="mediator")

summary(med.out)
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.01613      0.01872        0.013  <2e-16 *** 
# ADE             0.00214      0.00020        0.004   0.022 * 
# Total Effect    0.01827      0.02096        0.015  <2e-16 ***

results <- data.frame(
  EffectType = c("ACME", "ADE", "Total Effect"),
  Estimate = c(0.01613, 0.00214, 0.01827),
  LowerCI = c(0.01872, 0.00020, 0.02096),
  UpperCI = c(0.013, 0.004, 0.015)
)

pdf("./MM_FMT_RNAseq_mediation_test_D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae;D_5__uncultured bacterium.pdf", width = 8, height = 6)
ggplot(results, aes(x = EffectType, y = Estimate, fill = EffectType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(.9)) +
  theme_minimal() +
  labs(title = "Causal Mediation Analysis Results - D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae;D_5__uncultured bacterium",
       x = "Effect Type",
       y = "Effect Estimate") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank())
dev.off()

# Feature.ID: 850e0d58c1aa35238538564a8010e0ff

myData <- data.frame(mediator = mediator[,"850e0d58c1aa35238538564a8010e0ff"],
                     treatment = treatment,
                     outcome = outcome)

# Model for the mediator
mediator.model <- lm(mediator ~ treatment, data = myData)

# Model for the outcome
outcome.model <- lm(outcome ~ treatment + mediator, data = myData)

med.out <- mediate(mediator.model, outcome.model, treat="treatment", mediator="mediator")

summary(med.out)
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.016147     0.01872        0.013  <2e-16 ***  
# ADE             0.001946     0.00020        0.004    0.02 *  
# Total Effect    0.018093     0.02096        0.015  <2e-16 ***

results <- data.frame(
  EffectType = c("ACME", "ADE", "Total Effect"),
  Estimate = c(0.016147, 0.001946, 0.018093),
  LowerCI = c(0.01872, 0.00020, 0.02096),
  UpperCI = c(0.013, 0.004, 0.015)
)

pdf("./MM_FMT_RNAseq_mediation_test_D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae;D_5__uncultured bacterium;D_6__;D_7__.pdf", width = 8, height = 6)
ggplot(results, aes(x = EffectType, y = Estimate, fill = EffectType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(.9)) +
  theme_minimal() +
  labs(title = "Causal Mediation Analysis Results - D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae;D_5__uncultured bacterium;D_6__;D_7__",
       x = "Effect Type",
       y = "Effect Estimate") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank())
dev.off()

# Feature.ID: 55ebdcf041fb6b8b4632993423011372

myData <- data.frame(mediator = mediator[,"55ebdcf041fb6b8b4632993423011372"],
                     treatment = treatment,
                     outcome = outcome)

# Model for the mediator
mediator.model <- lm(mediator ~ treatment, data = myData)

# Model for the outcome
outcome.model <- lm(outcome ~ treatment + mediator, data = myData)

med.out <- mediate(mediator.model, outcome.model, treat="treatment", mediator="mediator")

summary(med.out)  
#                 Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.016123         0.02         0.01  <2e-16 ***
# ADE             0.002677        0.005        0.001   0.024 *  
# Total Effect    0.018800        0.024        0.015  <2e-16 ***

results <- data.frame(
  EffectType = c("ACME", "ADE", "Total Effect"),
  Estimate = c(0.016123, 0.002677, 0.018800),
  LowerCI = c(0.02, 0.005, 0.024),
  UpperCI = c(0.01, 0.001, 0.015))

pdf("./MM_FMT_RNAseq_mediation_test_D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae.pdf", width = 8, height = 6)
ggplot(results, aes(x = EffectType, y = Estimate, fill = EffectType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(.9)) +
  theme_minimal() +
  labs(title = "Causal Mediation Analysis Results - D_0__Bacteria;D_1__Bacteroidetes;D_2__Bacteroidia;D_3__Bacteroidales;D_4__Muribaculaceae",
       x = "Effect Type",
       y = "Effect Estimate") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank())
dev.off()

# Feature.ID: 54b3e0af28071e24c4f64e225e3b016d

myData <- data.frame(mediator = mediator[,"54b3e0af28071e24c4f64e225e3b016d"],
                     treatment = treatment,
                     outcome = outcome)

# Model for the mediator
mediator.model <- lm(mediator ~ treatment, data = myData)

# Model for the outcome
outcome.model <- lm(outcome ~ treatment + mediator, data = myData)

med.out <- mediate(mediator.model, outcome.model, treat="treatment", mediator="mediator")

summary(med.out) 
#                Estimate 95% CI Lower 95% CI Upper p-value    
# ACME            0.00269      0.00538         0.00   0.014 *  
# ADE             0.01344      0.01653         0.01  <2e-16 ***
# Total Effect    0.01613      0.01866         0.01  <2e-16 ***
# Prop. Mediated  0.16355      0.02202         0.32   0.014 *  

results <- data.frame(
  EffectType = c("ACME", "ADE", "Total Effect"),
  Estimate = c(0.00269, 0.01344, 0.01613),
  LowerCI = c(0.00538, 0.01653, 0.01866),
  UpperCI = c(0.00, 0.01, 0.01))

pdf("./MM_FMT_RNAseq_mediation_test_D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Ruminococcaceae;D_5__Ruminiclostridium 6;D_6__uncultured bacterium.pdf", width = 8, height = 6)
ggplot(results, aes(x = EffectType, y = Estimate, fill = EffectType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(.9)) +
  theme_minimal() +
  labs(title = "Causal Mediation Analysis Results - D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Ruminococcaceae;D_5__Ruminiclostridium 6;D_6__uncultured bacterium",
       x = "Effect Type",
       y = "Effect Estimate") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank())
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_mediation_analysis.txt", sep =""))
sessionInfo()
sink()
