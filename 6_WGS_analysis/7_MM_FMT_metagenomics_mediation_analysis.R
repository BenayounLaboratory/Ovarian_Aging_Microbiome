setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2025-01-20_Mediation")

library('mediation')
library(tidyr)
library(dplyr)
library(ggplot2)

######################################
# Mediation analysis:
#                     - Mediators: features (from FMT metagenomics data)
#                     - Treatment: FMT_YF vs. FMT_EF
#                     - Outcome: FMT ovarian bulk RNAseq
# Removed FMT_YF_7: High contamination %
# https://cran.r-project.org/web/packages/mediation/mediation.pdf
######################################

######################################
# 1. Import data
######################################

# Metagenomics data

FMT.Bracken <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/Kraken_outputs/MM_bracken_combined.out", 
                          header = TRUE, sep = "\t")
FMT.metadata <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2024-12-02/FMT_Metagenomics_metadata.txt", 
                           header = TRUE)
FMT.metadata$Group <- factor(FMT.metadata$Group, levels = c("YF", "OF"))

my.FMT.metagenomics.differentials <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2025-01-20/2025-01-20_MM_FMT_metagenomics_ALDEx2_differential_abundance_analysis_results.txt", 
                                                 header = TRUE, sep = "\t")

# Bulk ovarian RNA-seq

load("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2025-01-20_Mediation/2024-03-29_MM_FMT_Ovaries_DESeq2_Analysis_STAR_post_SVA_normalized_counts.RData")

tissue.cts.cl <- tissue.cts[,colnames(tissue.cts) != "FMT_YF_7"]

######################################
# 2. Pre-process data
######################################

# Filter significant features
my.FMT.metagenomics.differentials.sorted <- my.FMT.metagenomics.differentials[order(my.FMT.metagenomics.differentials$effect), ]
my.FMT.metagenomics.sig <- my.FMT.metagenomics.differentials.sorted[c(1:50),]

significant_features <- my.FMT.metagenomics.sig$Feature.ID

# Extract data from Bracken output

FMT.Bracken.sig <- FMT.Bracken[FMT.Bracken$name %in% significant_features, ]

# Filter only columns with `_num` in their names
FMT.Bracken.sig.num <- FMT.Bracken.sig[, grep("_num$", colnames(FMT.Bracken.sig))]

# Rename columns to remove everything after `MHK###`
colnames(FMT.Bracken.sig.num) <- gsub("\\.bracken_num$", "", colnames(FMT.Bracken.sig.num))
rownames(FMT.Bracken.sig.num) <- FMT.Bracken.sig$name

print(FMT.Bracken.sig.num)

# Perform CLR transformation
my.FMT.clr <- apply(log2(FMT.Bracken.sig.num+0.5), 2, function(x) x-mean(x))

# Treatment variable
treatment <- FMT.metadata$Group

# Mediator variable
t.my.FMT.table <- t(my.FMT.clr)

# Re-order rows to match
rownames(t.my.FMT.table) <- c("FMT_YF_1", "FMT_YF_2", 
                              "FMT_OF_1", "FMT_OF_2", "FMT_OF_3", 
                              "FMT_YF_3", "FMT_YF_4", "FMT_YF_5", "FMT_YF_6",
                              "FMT_OF_4", "FMT_OF_5", "FMT_OF_6", "FMT_OF_7")

rownames <- rownames(t.my.FMT.table)

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

t.my.FMT.table <- t.my.FMT.table[ordered_rownames, ]

######################################
# 3. Perform mediation analysis
######################################

mediator <- t.my.FMT.table

# Outcome variable - Use MDS 1 coordinates

mds.result <- cmdscale(1-cor(tissue.cts.cl, method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
outcome <- mds.result[, 1]

# Initialize an empty list to store significant results
significant_results <- list()

for (feature in colnames(mediator)) {

  myData <- data.frame(mediator = mediator[, feature],
                       treatment = treatment,
                       outcome = outcome)
  
  # Model for mediator
  mediator.model <- lm(mediator ~ treatment, data = myData)
  
  # Model for outcome
  outcome.model <- lm(outcome ~ treatment + mediator, data = myData)
  
  # Perform mediation analysis
  med.out <- mediate(mediator.model, outcome.model, treat = "treatment", mediator = "mediator")
  
  # Summarize the mediation analysis results
  summary <- summary(med.out)
  
  # Check if ACME p-value is significant
  if (summary$d0.p < 0.05) {
    # Store the results in a data frame
    results <- data.frame(
      EffectType = c("ACME", "ADE", "Total Effect"),
      Estimate = c(summary$d0, summary$z0, summary$tau.coef),
      LowerCI = c(summary$d0.ci[1], summary$z0.ci[1], summary$tau.ci[1]),
      UpperCI = c(summary$d0.ci[2], summary$z0.ci[2], summary$tau.ci[2]),
      PValue = c(summary$d0.p, summary$z0.p, summary$tau.p),
      Feature = feature
    )
    
    # Save the results into the list of significant results
    significant_results[[feature]] <- results
    
    # Plot results for the current significant feature
    pdf(paste0("./MM_FMT_RNAseq_mediation_test_", feature, ".pdf"), width = 8, height = 6)
    ggplot(results, aes(x = EffectType, y = Estimate, fill = EffectType)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(.9)) +
      theme_minimal() +
      labs(
        title = paste("Causal Mediation Analysis Results -", feature),
        x = "Effect Type",
        y = "Effect Estimate"
      ) +
      scale_fill_brewer(palette = "Set1") +
      theme(legend.title = element_blank())
    dev.off()
  }
}

save(significant_results, file = paste0(Sys.Date(), "_MM_FMT_metagenomics_mediation_analysis_results_stats_summary.RData"))

###########################################################
sink(file = paste(Sys.Date(),"MM_metagenomics_mediation_analysis_session_info.txt", sep =""))
sessionInfo()
sink()
