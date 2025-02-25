setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Ovarian_health_index/")

library(dplyr)
library('beeswarm')
library(outliers)

###################################
# Menopause-microbiome project - revision
# Calculate ovarian health index scores
# Parameters: AMH, FSH, INHBA
#             Follicle count (combined)
# Cohorts: AC females, VCD, FMT
###################################

####################################
# 1. Import data files
###################################

# AC - females

my.hormone.data.AC.females <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Input_files/MM_AC_cohort_females_hormone_data.txt", header= TRUE, sep = "\t")
my.follicle.data.AC.females <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Input_files/MM_AC_cohort_females_ovarian_follicle_count_median.txt", header= TRUE, sep = "\t")

my.follicle.data.AC.females$combined <- rowSums(my.follicle.data.AC.females[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

# VCD cohort

my.hormone.data.VCD <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Input_files/MM_VCD_cohort_hormone_data.txt", header= TRUE, sep = "\t")
my.follicle.data.VCD <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Input_files/MM_VCD_cohort_ovarian_follicle_count_median.txt", header= TRUE, sep = "\t")

my.follicle.data.VCD$combined <- rowSums(my.follicle.data.VCD[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

# FMT cohort

my.hormone.data.FMT <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Serum_hormone_quantification_FMT/2025-02-10_MeMo_hormone_data_FMT_corrected.txt", header= TRUE, sep = "\t")
my.follicle.data.FMT <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Follicle_counts_FMT/MM_FMT_follicle_counts_combined.txt", header= TRUE, sep = "\t")

my.follicle.data.FMT$combined <- rowSums(my.follicle.data.FMT[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

####################################
# 2. Assign scores based on AC data
###################################

# AC data

my.AC.data <- cbind(my.hormone.data.AC.females, my.follicle.data.AC.females[, c("combined")])
colnames(my.AC.data)[ncol(my.AC.data)] <- "follicle_combined"

# VCD cohort

my.VCD.data <- cbind(my.hormone.data.VCD, my.follicle.data.VCD[, c("combined")])
colnames(my.VCD.data)[ncol(my.VCD.data)] <- "follicle_combined"

# FMT cohort

my.FMT.data <- cbind(my.hormone.data.FMT, my.follicle.data.FMT[, c("combined")])
colnames(my.FMT.data)[ncol(my.FMT.data)] <- "follicle_combined"

colnames(my.FMT.data) <- c("Sample_ID", "Mouse_ID", "Treatment", "AMH", "FSH_Multiplex", "FSH_US", "INHBA", "Cohort", "FSH_normalized", "FSH", "follicle_combined")

# Calculate medians for YF and EF groups
medians <- my.AC.data %>%
  group_by(Group) %>%
  summarise(across(c(AMH, FSH, INHBA, follicle_combined), median, .names = "{.col}_median"))

YF_medians <- filter(medians, Group == "YF")
EF_medians <- filter(medians, Group == "EF")

# Function to assign scores
calculate_scores <- function(value, YF_median, EF_median, increasing) {
  if (increasing) {
    return(case_when(
      value <= YF_median ~ 3,
      value > YF_median & value < EF_median ~ 2,
      TRUE ~ 1
    ))
  } else {
    return(case_when(
      value <= EF_median ~ 1,
      value > EF_median & value < YF_median ~ 2,
      TRUE ~ 3
    ))
  }
}

# Function to process datasets and compute scores
process_data <- function(data) {
  data <- data %>%
    rowwise() %>%
    mutate(
      AMH_score = calculate_scores(AMH, YF_medians$AMH_median, EF_medians$AMH_median, FALSE),
      FSH_score = calculate_scores(FSH, YF_medians$FSH_median, EF_medians$FSH_median, TRUE),
      INHBA_score = calculate_scores(INHBA, YF_medians$INHBA_median, EF_medians$INHBA_median, FALSE),
      hormone_avg_score = (AMH_score + FSH_score + INHBA_score) / 3,
      follicle_score = calculate_scores(follicle_combined, YF_medians$follicle_combined_median, EF_medians$follicle_combined_median, FALSE),
      total_score = (hormone_avg_score + follicle_score) / 6 * 100
    )
  return(data)
}

# Process all datasets
my.AC.data <- process_data(my.AC.data)
my.VCD.data <- process_data(my.VCD.data)
my.FMT.data <- process_data(my.FMT.data)

# Save processed data
write.table(my.AC.data, "MM_Ovarian_health_index_AC_cohort.txt", quote = FALSE, sep = "\t")
write.table(my.VCD.data, "MM_Ovarian_health_index_VCD_cohort.txt", quote = FALSE, sep = "\t")
write.table(my.FMT.data, "MM_Ovarian_health_index_FMT_cohort.txt", quote = FALSE, sep = "\t")

####################################
# 3. Generate boxplots
###################################

generate_boxplot <- function(data, group_col, groups, colors, title, filename) {
  index_list <- setNames(
    list(
      data$total_score[data[[group_col]] == groups[1]],
      data$total_score[data[[group_col]] == groups[2]]
    ),
    groups
  )
  
  p_value <- wilcox.test(index_list[[1]], index_list[[2]])$p.value
  
  pdf(filename, width = 4, height = 7)
  boxplot(index_list, 
          outline = FALSE, ylim = c(0, 120), 
          col = colors, las = 1, 
          ylab = "Ovarian health index", main = title)
  beeswarm(index_list, pch = 16, col = "black", add = TRUE, cex = 1)
  text(1.5, 10, paste("p-value =", format(p_value, digits = 5)))
  dev.off()
}

generate_boxplot(my.AC.data, "Group", c("YF", "EF"), c("deeppink1", "deeppink4"),
                 "AC females - Ovarian health index", "MM_AC_cohort_females_ovarian_index_score_boxplot.pdf")

generate_boxplot(my.VCD.data, "Treatment", c("CTL", "VCD"), c("green", "yellow"),
                 "VCD cohort - Ovarian health index", "MM_VCD_cohort_ovarian_index_score_boxplot.pdf")

generate_boxplot(my.FMT.data, "Treatment", c("FMT_YF", "FMT_EF"), c("purple", "orange"),
                 "FMT cohort - Ovarian health index", "MM_FMT_cohort_ovarian_index_score_boxplot.pdf")

####################################
# 4. Replot FMT cohorts - distinguish two cohorts (different FSH kits used)
###################################

my.FMT.data$Treatment <- factor(my.FMT.data$Treatment, levels = c("FMT_YF","FMT_EF"))

# Define the dataset for Wilcoxon test
FMT.data.list <- list(
  "FMT_YF" = my.FMT.data$total_score[my.FMT.data$Treatment == "FMT_YF"],
  "FMT_EF" = my.FMT.data$total_score[my.FMT.data$Treatment == "FMT_EF"]
)

# Perform Wilcoxon test
p_value <- wilcox.test(FMT.data.list[["FMT_YF"]], FMT.data.list[["FMT_EF"]])$p.value

# Extract unique cohort levels and assign different pch values
cohort_levels <- unique(my.FMT.data$Cohort)
pch_values <- 1:length(cohort_levels) 

# Create a mapping of cohort to pch values
cohort_pch_map <- setNames(pch_values, cohort_levels)

# Assign pch values based on cohort
my.FMT.data$pch_values <- cohort_pch_map[my.FMT.data$Cohort]

# Define colors for the groups
box_colors <- c("purple", "orange")

# Generate the boxplot
pdf(paste0(Sys.Date(), "_MM_FMT_ovarian_health_index_by_cohort.pdf"), width = 4, height = 7)

boxplot(FMT.data.list,
        outline = FALSE, ylim = c(0, 120),  
        col = box_colors, las = 1,
        ylab = "Ovarian Health Index", main = "FMT Effect on Ovarian Health")
beeswarm(total_score ~ Treatment, 
         data = my.FMT.data,
         col = "black", 
         pwpch = my.FMT.data$pch_values, 
         add = TRUE, cex = 1.2) 
legend("topright", legend = cohort_levels, pch = pch_values, col = "black", title = "Cohort")
mtext(paste("Wilcoxon p-value:", signif(p_value, 3)), side = 3, line = -1.5, cex = 0.9)

dev.off()

################################################################################
sink(file = paste(Sys.Date(), "_MM_revision_calculate_ovarian_health_index_updated_Session_Info.txt"))
sessionInfo()
sink()
