setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/")

library(mia) 
library(SummarizedExperiment)
library(phyloseq)
library(vegan)
library(picante)
library(ecodist)
library(ggplot2)
library(compositions)
library(mixOmics)
library(dplyr)
library(beeswarm)
library(ALDEx2)
library(tidyverse)
library(ggrepel)


rm(list=ls())

#############################################################
# MM revision - FMT metagenomics
# Analyze Bracken outfiles - batch correction & PCA analysis
# Remove MHK271 -> high host DNA alignment rate (contamination)
#############################################################

#############################################################
# 1. Load data
#############################################################

FMT.Bracken <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/Kraken_outputs/MM_bracken_combined.out", 
                          header = TRUE, sep = "\t")
FMT.metadata <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2024-12-02/FMT_Metagenomics_metadata.txt", 
                           header = TRUE)

# Remove MHK271
FMT.Bracken.cl <- FMT.Bracken[, !grepl("MHK271", colnames(FMT.Bracken))]
FMT.metadata.cl <- FMT.metadata[FMT.metadata$Sample_ID != c("MHK271", "MHK272"),]

#############################################################
# 2. QC and Filtering
#############################################################

# Identify count and fraction columns
count_cols <- grep("bracken_num", colnames(FMT.Bracken.cl))
frac_cols <- grep("bracken_frac", colnames(FMT.Bracken.cl))

# Filter by count-based threshold
count_threshold <- 20
FMT.Bracken.count.filtered <- FMT.Bracken.cl[rowSums(FMT.Bracken.cl[, count_cols]) > count_threshold, ]

# Filter by relative abundance 
frac_threshold <- 0.01 / 100 
FMT.Bracken.frac.filtered <- FMT.Bracken.count.filtered[rowMeans(FMT.Bracken.count.filtered[, frac_cols]) > frac_threshold, ]

# Filter by prevalence 
prevalence_threshold <- 0.10 * length(count_cols) 
FMT.Bracken.final.filtered <- FMT.Bracken.frac.filtered[rowSums(FMT.Bracken.frac.filtered[, count_cols] > 0) >= prevalence_threshold, ]

FMT.Bracken.filtered <- FMT.Bracken.final.filtered     # 273 species detected
rownames(FMT.Bracken.filtered) <- FMT.Bracken.filtered$name

count_data <- FMT.Bracken.filtered[,count_cols]
frac_data   <- FMT.Bracken.filtered[,frac_cols]

# Cleanup column names
colnames(count_data) <- gsub("\\..*", "", colnames(count_data))
colnames(frac_data) <- gsub("\\..*", "", colnames(frac_data))

colnames(count_data)
colnames(frac_data)

tax_ids <- as.matrix(FMT.Bracken.filtered$taxonomy_id)
write.table(tax_ids, file = paste0(Sys.Date(), "_MM_FMT_Metagenomics_tax_IDs.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#############################################################
# 2. PCA - CLR transformed data
#############################################################

# CLR transformation on filtered data
filtered_data <- frac_data + 1e-6 
clr_data_filtered <- clr(frac_data)
clr_data_filtered <- as.data.frame(clr_data_filtered)

# Perform PCA

pca_result <- pca(t(clr_data_filtered), ncomp = 3)
pca_result$variates$X

# Add Group metadata as needed
clr_data_filtered_with_metadata <- rbind(clr_data_filtered, Group = FMT.metadata$Group[match(colnames(clr_data_filtered), FMT.metadata$Sample_ID)])

# Perform PCA
pca_result <- pca(t(clr_data_filtered), ncomp = 3)

# Extract the scores for the first three principal components
pca_data <- as.data.frame(pca_result$x[, 1:3])
colnames(pca_data) <- c("PC1", "PC2", "PC3")

# Add metadata for visualization
# Extract sample IDs from row names in pca_data
pca_data$Sample_ID <- sub("\\..*", "", rownames(pca_data))
pca_data_with_group <- merge(pca_data, FMT.metadata, by = "Sample_ID", all.x = TRUE)

pca_data_with_group <- pca_data_with_group %>% select(Sample_ID, Group, PC1, PC2, PC3)

# Plot PCA results
pdf(paste0(Sys.Date(), "_CLR_PCA_plot.pdf"), width = 7, height = 7)
ggplot(pca_data_with_group, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of CLR-transformed Data") +
  theme_minimal()
dev.off()

write.table(clr_data_filtered_with_metadata, paste0(Sys.Date(), "_MM_FMT_metagenomics_CLR_transformed_data.txt"), sep = "\t", quote = FALSE)

#############################################################
# 3. Beta diversity analysis - Bray-Curtis & Jaccard
#############################################################

# Function to calculate PCoA and plot results
perform_pcoa_plot <- function(distance_method, title) {
  
  dist_matrix <- vegan::vegdist(t(frac_data), method = distance_method)
  
  pcoa_result <- ecodist::pco(dist_matrix)
  
  pcoa_df <- data.frame(
    pcoa1 = pcoa_result$vectors[,1],
    pcoa2 = pcoa_result$vectors[,2],
    Sample_ID = sub("\\.bracken_frac$", "", rownames(pcoa_result$vectors))
  )
  
  pcoa_df <- merge(pcoa_df, FMT.metadata, by = "Sample_ID")
  
  plot <- ggplot(data = pcoa_df, aes(x = pcoa1, y = pcoa2, color = Group)) +
    geom_point(size = 3) +
    labs(x = "PC1",
         y = "PC2", 
         title = title) +
    theme(title = element_text(size = 10))
  
  return(plot)
}

bray_curtis_plot <- perform_pcoa_plot("bray", "Bray-Curtis PCoA")
jaccard_plot <- perform_pcoa_plot("jaccard", "Jaccard PCoA")

pdf(paste0(Sys.Date(), "_Bray_curtis_plot.pdf"), width = 7, height = 7)
print(bray_curtis_plot)
dev.off()

pdf(paste0(Sys.Date(), "_Jaccard_plot.pdf"), width = 7, height = 7)
print(jaccard_plot)
dev.off()

#############################################################
# 4. Measure alpha diversity
#############################################################

FMT.metadata.cl <- as.data.frame(FMT.metadata.cl)

phy_tree <- read.tree("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/FMT_metagenomics/2024-12-02/Phylo_tree/phyliptree.phy")

# Clean tree tip labels
phy_tree$tip.label <- gsub("^'|'$", "", gsub(" ", "_", phy_tree$tip.label))

count_matrix <- as.matrix(count_data)
rownames(count_matrix) <- gsub(" ", "_", rownames(count_matrix))

# Ensure compatibility between tree and count matrix
stopifnot(all(rownames(count_matrix) %in% phy_tree$tip.label)) 

# Create the TreeSummarizedExperiment object
tse <- TreeSummarizedExperiment(
  assays = list(counts = count_matrix),
  rowTree = phy_tree 
)

# Clean and add metadata
colnames(tse) <- gsub("\\.bracken_num$", "", colnames(tse)) 
rownames(FMT.metadata.cl) <- FMT.metadata.cl$Sample_ID      
FMT.metadata.cl <- FMT.metadata.cl[, -1]                 
colData(tse) <- DataFrame(FMT.metadata.cl)   

# Calculate diversity metrics
tse <- estimateDiversity(tse, index = "shannon", assay.type = "counts")
tse <- estimateDiversity(tse, index = "faith", assay.type = "counts", tree = "phylo") 

# Calculate observed features
colData(tse)$observed <- colSums(assay(tse, "counts") > 0) 

# View updated diversity metrics in colData
colData(tse)

df <- as.data.frame(colData(tse))

observed.features <- list(
  "YF" = df$observed[df$FMT.metadata.cl == "YF"],
  "OF" = df$observed[df$FMT.metadata.cl == "OF"]
)

shannon.entropy <- list(
  "YF" = df$shannon[df$FMT.metadata.cl == "YF"],
  "OF" = df$shannon[df$FMT.metadata.cl == "OF"]
)

# Perform Wilcoxon test 
observed_wilcox <- wilcox.test(df$observed[df$FMT.metadata.cl == "YF"], 
                               df$observed[df$FMT.metadata.cl == "OF"])
shannon_wilcox <- wilcox.test(df$shannon[df$FMT.metadata.cl == "YF"], 
                              df$shannon[df$FMT.metadata.cl == "OF"])

pdf(paste0(Sys.Date(), "_Observed_features.pdf"), width = 4, height = 7)
boxplot(observed.features, 
        outline = TRUE, 
        col = c("purple", "orange"),
        main = "Observed Features", 
        ylab = "Observed Features", 
        xlab = "Group",
        ylim = c(min(df$observed) - 5, max(df$observed) + 5))
beeswarm(observed.features, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5, max(df$observed) + 3, paste("p-value =", signif(observed_wilcox$p.value, 3)), cex = 1)
dev.off()

pdf(paste0(Sys.Date(), "_Shannon_entropy.pdf"), width = 4, height = 7)
boxplot(shannon.entropy, 
        outline = TRUE, 
        col = c("purple", "orange"),
        main = "Shannon Entropy", 
        ylab = "Shannon Entropy", 
        xlab = "Group",
        ylim = c(min(df$shannon) - 0.5, max(df$shannon) + 0.5))
beeswarm(shannon.entropy, pch = 16, col = "black", add = TRUE, cex = 1.0)
text(1.5, max(df$shannon) + 0.3, paste("p-value =", signif(shannon_wilcox$p.value, 3)), cex = 1)
dev.off()

#############################################################
# 5. Differential abundance analysis - use CLR-transformed data
#############################################################

clr_data_filtered <- as.data.frame(clr_data_filtered)

FMT.metadata.cl <- FMT.metadata[FMT.metadata$Sample_ID %in% colnames(count_data), ]
group <- FMT.metadata.cl$Group  
names(group) <- FMT.metadata.cl$Sample_ID

group <- as.character(group)
count_data <- as.matrix(count_data)

# ALDEx2 analysis
aldex_res <- aldex.clr((count_data), group, mc.samples = 128, denom = "all", verbose = TRUE)
aldex_effect <- aldex.ttest(aldex_res)
aldex_summary <- aldex.effect(aldex_res)

aldex_combined <- data.frame(
  Feature.ID = rownames(aldex_effect),
  aldex_effect,
  aldex_summary
)

# Filter significant features
aldex_combined$significant <- ifelse(aldex_combined$wi.ep < 0.05, "Significant", "Not Significant")

# Replace 0 values with a small number
aldex_combined$wi.ep[aldex_combined$wi.ep == 0] <- 1e-300

# Cap the maximum -log10(p-value)
cap_value <- 3 
aldex_combined$minus_log10_pvalue <- pmin(-log10(aldex_combined$wi.ep), cap_value)

# Replace 0 p-values for log transformation
aldex_combined$wi.ep[aldex_combined$wi.ep == 0] <- 1e-300

# Calculate -log10(pval)
cap_value <- 3 
aldex_combined$minus_log10_pvalue <- pmin(-log10(aldex_combined$wi.ep), cap_value)

aldex_combined$color <- ifelse(
  aldex_combined$wi.ep < 0.05,  
  aldex_combined$effect,     
  NA       
)

# Volcano Plot
aldex_combined <- aldex_combined %>%
  mutate(
    label = ifelse(
      Feature.ID %in% c(
        rownames(aldex_combined)[order(effect, decreasing = TRUE)[1:5]],
        rownames(aldex_combined)[order(effect, decreasing = FALSE)[1:5]]
      ),
      Feature.ID,
      NA
    )
  )

pdf(paste0(Sys.Date(), "_FMT_metagenomics_differential_abundance_volcano_plot_labeled.pdf"), width = 8, height = 4)
ggplot(aldex_combined, aes(x = effect, y = minus_log10_pvalue)) +
  geom_point(aes(color = color), alpha = 0.8, size = 3) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "gray" 
  ) +
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = 10
  ) +
  labs(
    title = "Volcano Plot of Differential Features",
    x = "Effect Size",
    y = "-log10(P-value)",
    color = "Effect Size (Significant Only)"
  ) +
  theme_minimal()
dev.off()

# Bubble plot for top 10 and bottom 10 features based on effect size
top_bottom_features <- aldex_combined %>%
  arrange(desc(effect)) %>%
  slice_head(n = 10) %>%  
  bind_rows(
    aldex_combined %>%
      arrange(effect) %>%
      slice_head(n = 10)
  )

pdf(paste0(Sys.Date(), "_Bubble_Plot_Top_Bottom_10_Features.pdf"), width = 8, height = 6)
ggplot(top_bottom_features, aes(x = reorder(Feature.ID, effect), y = effect)) +
  geom_point(aes(size = minus_log10_pvalue, color = effect), alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "gray" 
  ) +
  scale_size_continuous(range = c(3, 10)) +  
  coord_flip() +  
  labs(
    title = "Bubble Plot of Top and Bottom 10 Features",
    x = "Feature ID",
    y = "Effect Size",
    size = "-log10(P-value)",
    color = "Effect Size"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
dev.off()

# Extract the list of top and bottom 10 feature IDs
top_bottom_feature_ids <- top_bottom_features$Feature.ID

selected_count_data <- count_data[top_bottom_feature_ids, , drop = FALSE]

long_count_data <- as.data.frame(t(selected_count_data))
long_count_data$Sample_ID <- rownames(long_count_data) 

long_count_data <- reshape2::melt(
  long_count_data,
  id.vars = "Sample_ID",
  variable.name = "Feature_ID",
  value.name = "Count"
)

# Add metadata 
long_count_data <- merge(long_count_data, FMT.metadata, by = "Sample_ID", all.x = TRUE)

long_count_data$Feature_ID <- factor(
  long_count_data$Feature_ID,
  levels = top_bottom_feature_ids
)

unique_features <- levels(long_count_data$Feature_ID)

# Create boxplots 
pdf(paste0(Sys.Date(), "_MM_FMT_cohort_base_abundance_boxplots.pdf"), width = 15, height = 15)

num_features <- length(unique_features) 
par(mfrow = c(ceiling(num_features / 4), 4), mar = c(4, 4, 2, 1))  

for (feature in unique_features) {

  feature_data <- subset(long_count_data, Feature_ID == feature)
  
  feature_data$Group <- factor(feature_data$Group, levels = c("YF", "OF"))
  
  boxplot(
    Count ~ Group,
    data = feature_data,
    main = as.character(feature), 
    xlab = "FMT Group",
    ylab = "Counts",
    col = c("purple", "orange"), 
    border = "black", 
    outline = TRUE 
  )
  
  beeswarm::beeswarm(
    Count ~ Group,
    data = feature_data,
    add = TRUE,  
    pch = 16, 
    col = "black",  
    method = "center",  
    cex = 0.8 
  )
}

par(mfrow = c(1, 1))
dev.off()

write.table(aldex_combined, paste0(Sys.Date(), "_MM_FMT_metagenomics_ALDEx2_differential_abundance_analysis_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(unique_features, paste0(Sys.Date(), "_MM_FMT_metagenomics_sig_features.txt"), quote = FALSE, sep = "\t")

################################################################################
sink(file = paste(Sys.Date(), "_MM_FMT_metagenomics_analysis_Session_Info.txt"))
sessionInfo()
sink()