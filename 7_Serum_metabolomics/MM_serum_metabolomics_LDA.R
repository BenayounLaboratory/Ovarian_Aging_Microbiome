setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Analysis_V6")

library(MASS) 
library(tidyverse)
library(beeswarm)

############################################
# Menopause-microbiome project - revision
# Serum metabolomics - LDA anlaysis
############################################

############################################
# 1. Load data & pre-process data
############################################

################
# Load matrix & ID mapping data

data <- read.csv('/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Raw_data/fc_log_fdr_peak_ms2_filtered_matrix.csv', row.names = 1)
ID.mapping <- read.csv('/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Raw_data/ID_to_chemical_name_mapping.csv', row.names = 1)

################
# Load metadata and combine with matrix

metadata <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Metabolomics_serum/Raw_data/MM_metabolomics_serum_metadata.txt", header = TRUE, sep = "\t")

# Merge metadata with metabolomics data
data$Sample_ID <- rownames(data)
merged_data <- left_join(metadata, data, by = "Sample_ID")

# Convert Group and Batch to factors
merged_data$Group <- as.factor(merged_data$Group)
merged_data$Batch <- as.factor(merged_data$Batch)

# Identify contradictory metabolites

# Compute log2 fold change within each batch
fc_data <- merged_data %>%
  pivot_longer(cols = -c(Sample_ID, Group, Batch), names_to = "Metabolite", values_to = "Value") %>%
  group_by(Batch, Metabolite, Group) %>%
  summarise(Median_Value = median(Value, na.rm = TRUE), .groups = "drop") %>%  # Compute median per batch and group
  pivot_wider(names_from = Group, values_from = Median_Value) %>%  # Convert FMT_YF and FMT_EF into columns
  mutate(log2FC = log2(FMT_EF / FMT_YF))  # Compute log2 fold change

# Identify metabolites with opposite trends between batches
fc_wide <- fc_data %>%
  select(Batch, Metabolite, log2FC) %>%
  pivot_wider(names_from = Batch, values_from = log2FC, names_prefix = "Batch")

# Find metabolites where FC directions are opposite across batches
opposing_metabolites <- fc_wide %>%
  filter((Batch1 > 0 & Batch2 < 0) | (Batch1 < 0 & Batch2 > 0)) %>%
  pull(Metabolite)

# print(opposing_metabolites)

############################################
# 2. Perform LDA
############################################

# Scale data

data <- as.data.frame(scale(data))

# Merge metadata with metabolomics data
data$Sample_ID <- rownames(data)
merged_data <- left_join(metadata, data, by = "Sample_ID")

# Convert Group and Batch to factors
merged_data$Group <- as.factor(merged_data$Group)
merged_data$Batch <- as.factor(merged_data$Batch)

# Remove non-numeric columns

data.clean <- merged_data %>% select(-Sample_ID, -Batch)

# Perform LDA
lda_model <- lda(Group ~ ., data = data.clean)

# Print model summary
print(lda_model)

# Extract LDA scores
lda_scores <- predict(lda_model)$x

# Plot LDA distribution
lda_df <- data.frame(Sample_ID = merged_data$Sample_ID, Group = merged_data$Group, LDA1 = lda_scores[,1])

pdf(paste0(Sys.Date(), "_MM_serum_metabolomics_LDA_analysis_LDA_projection_scatter_plot.pdf"))
ggplot(lda_df, aes(x = LDA1, y = 0, color = Group)) +
  geom_point(size = 3, alpha = 0.8, position = position_nudge(y = 0)) + 
  scale_color_manual(values = c("FMT_YF" = "purple", "FMT_EF" = "orange")) + 
  labs(title = "LDA Projection of Serum Metabolomics Data",
       x = "LDA 1 Score", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(), 
        panel.grid.major.y = element_blank()) 
dev.off()


############################################
# 3. Assess LDA loadings
############################################

# Remove confounded metabolites
filtered_lda_loadings <- as.data.frame(lda_model$scaling[setdiff(rownames(lda_model$scaling), opposing_metabolites), , drop = FALSE])

# Sort by absolute value of LD1 
sorted_lda_loadings <- filtered_lda_loadings %>%
  as.data.frame() %>%
  mutate(Abs_LD1 = abs(LD1)) %>%
  arrange(desc(Abs_LD1)) %>%
  select(-Abs_LD1)

top_20_metabolites <- head(sorted_lda_loadings, 15)

# Retrieve original values 
metabolite_ids <- gsub("^X", "", rownames(top_20_metabolites))

# Map metabolite IDs to chemical names
top_20_metabolites$Chemical_Name <- ID.mapping[metabolite_ids, "Name"]
 
top_20_metabolites$Chemical_Name[is.na(top_20_metabolites$Chemical_Name)] <- metabolite_ids

# Plot top 20 metabolites with actual values
pdf(paste0(Sys.Date(), "_MM_serum_metabolomics_LDA_analysis_LDA_score_top_20_barplot.pdf"), height = 10, width = 15)
ggplot(top_20_metabolites, aes(x = reorder(Chemical_Name, -LD1), y = LD1, fill = LD1 > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("orange", "purple"), labels = c("Higher in FMT_EF", "Higher in FMT_YF")) +
  labs(title = "Top 20 Metabolites Contributing to LDA Separation",
       x = "Metabolites", y = "LDA Coefficient", fill = "Effect Direction") +
  theme_minimal()
dev.off()

# Subset the original data to include selected metabolites and Batch info
boxplot_data <- merged_data %>%
  select(Group, Batch, all_of(rownames(top_20_metabolites))) %>%
  pivot_longer(cols = -c(Group, Batch), names_to = "Metabolite", values_to = "Value")

boxplot_data$Group <- factor(boxplot_data$Group, levels = c("FMT_YF", "FMT_EF"))
boxplot_data$Metabolite <- factor(boxplot_data$Metabolite, levels = rownames(top_20_metabolites))

# Generate boxplots 

generate_metabolite_boxplot <- function(metabolite_name, y_limit) {
  # Subset the data for the metabolite
  df_metab <- boxplot_data %>%
    filter(Metabolite == metabolite_name)
  
  df_metab$Group <- factor(df_metab$Group, levels = c("FMT_YF", "FMT_EF"))
  
  # Perform Wilcoxon test
  wilcox_result <- wilcox.test(Value ~ Group, data = df_metab)
  print(paste("Wilcoxon p-value for", metabolite_name, ":", wilcox_result$p.value))
  
  pdf(paste0(Sys.Date(), "_", gsub(" ", "_", metabolite_name), "_metabolite_boxplot.pdf"), width = 5, height = 7)
  boxplot(Value ~ Group, data = df_metab,
          outline = FALSE,
          ylim = y_limit,
          col = c("purple", "orange"),
          las = 1,
          ylab = "Metabolite Abundance",
          main = metabolite_name)
  beeswarm(Value ~ Group, data = df_metab %>% filter(Batch == 1), 
           pch = 16, col = "black", add = TRUE, cex = 1.2, 
           method = "center", offset = -2) 
  
  beeswarm(Value ~ Group, data = df_metab %>% filter(Batch == 2), 
           pch = 17, col = "black", add = TRUE, cex = 1.2, 
           method = "center", offset = 2) 
    text(1.5, y_limit[2] * 0.95, paste("p-value ~", round(wilcox_result$p.value, 5)), cex = 0.8)
  dev.off()
}

for (metab in rownames(top_20_metabolites)) {
  generate_metabolite_boxplot(metab, y_limit = c(min(boxplot_data$Value, na.rm = TRUE), max(boxplot_data$Value, na.rm = TRUE)))
}

############################################
sink(file = paste0(Sys.Date(), "_MM_serum_metabolomics_LDA_analysis_session_info.txt"))
sessionInfo()
sink()
