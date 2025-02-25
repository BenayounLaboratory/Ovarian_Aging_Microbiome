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
# FMT cohort
###################################

###################################
# 1. Import picrust2 output files
###################################

# FMT cohort

FMT.pathway.pred <- read_qza("./q2-picrust2_output_FMT/pathway_abundance.qza")$data
FMT.metadata <- read_q2metadata("../../1_Manifest_and_metadata/MM_FMT_AC_PostFMT_sample-metadata.tsv")

# map file

map.file <- read.table("./metacyc_pathways_info.txt", sep = "\t")

# Convert values into integers

FMT.pathway.pred <- round(FMT.pathway.pred)

###################################
# 2. Perform ALDEx2 corr
###################################

set.seed(12345)

conds <- ifelse(grepl("Y", FMT.metadata$fmt), "FMT_YF", "FMT_EF")
conds.int <- ifelse(grepl("FMT_YF", conds), 1,2)

# Run Aldex2

aldex2 = aldex(FMT.pathway.pred, conds, mc.samples = 500, test = "t", 
               effect = TRUE, denom = "iqlr", verbose = TRUE)

head(aldex2, 10)

x <- aldex.clr(FMT.pathway.pred, conds)
corr.test <- aldex.corr(x, conds.int)

corr.test$log_p <- -log(corr.test$kendall.ep)
corr.test <- corr.test[order(-corr.test$log_p),]

# Convert row names of aldex2 to a column.
aldex2 <- aldex2 %>% 
  rownames_to_column(var = "Pathway")

# Perform a left join to map descriptions based on the Pathway column.
aldex2 <- aldex2 %>% 
  left_join(map.file, by = c("Pathway" = "V1"))

# Optionally, rename the V2 column to "description".
colnames(aldex2)[colnames(aldex2) == "V2"] <- "description"

aldex2_filtered <- subset(aldex2, effect >= 0.5)
aldex2_filtered <- aldex2_filtered[!is.na(aldex2_filtered$description), ]

aldex2_filtered$adjP <- -log10(aldex2_filtered$wi.ep)

aldex2_filtered_sorted <- aldex2_filtered[order(-aldex2_filtered$adjP), ]

pdf(paste0(Sys.Date(), "_MM_FMT_picrust2_barplot_down_in_FMT.pdf"), width = 15, height = 5)
ggplot(aldex2_filtered_sorted, aes(x = reorder(description, adjP), y = adjP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + # Flip coordinates for horizontal bars
  labs(x = "description", y = "-log10(adjP)", title = "Pathway vs. Adjusted P-value") +
  theme_minimal()
dev.off()

aldex2_filtered <- subset(aldex2, effect <= -0.5)
aldex2_filtered <- aldex2_filtered[!is.na(aldex2_filtered$description), ]

aldex2_filtered$adjP <- -log10(aldex2_filtered$wi.ep)

aldex2_filtered_sorted <- aldex2_filtered[order(-aldex2_filtered$adjP), ]

pdf(paste0(Sys.Date(), "_MM_FMT_picrust2_barplot_up_in_FMT.pdf"), width = 15, height = 5)
ggplot(aldex2_filtered_sorted, aes(x = reorder(description, adjP), y = adjP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + # Flip coordinates for horizontal bars
  labs(x = "description", y = "-log10(adjP)", title = "Pathway vs. Adjusted P-value") +
  theme_minimal()
dev.off()

###########################################################
sink(file = paste(Sys.Date(),"MM_FMT_cohort_functional_abundance_differential_abundance_analysis.txt", sep =""))
sessionInfo()
sink()
