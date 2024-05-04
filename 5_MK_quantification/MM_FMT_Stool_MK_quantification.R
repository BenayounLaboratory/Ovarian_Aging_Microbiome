library(ggplot2)
library(tidyverse)
library(beeswarm)

################################################################################
# Menopause_microbiome project
# MK quantification from Booth lab - Stool
################################################################################

vk2.quant.stool <- read.table("./2023-07-28_VK_quantification_Booth_lab_rawdata_stool.txt", sep = "\t", header = TRUE)
vk2.quant.stool <- vk2.quant.stool[c(1, 4:13)]

cols_to_sum <- colnames(vk2.quant.stool)
cols_to_sum <- cols_to_sum[2:10]

vk2.quant.stool$combined <- rowSums(vk2.quant.stool[, cols_to_sum])

# print(vk2.quant.stool)

###########################
# Stool - raw value
###########################

vk2.quant.stool.long <- vk2.quant.stool %>%
  pivot_longer(cols = MK5:MK13, names_to = "MK", values_to = "Value")

vk2.quant.stool.long$FMT <- factor(vk2.quant.stool.long$FMT, levels = c("FMT_YF", "FMT_OF"))
vk2.quant.stool.long$MK <- factor(vk2.quant.stool.long$MK, levels = c("MK5", "MK6", "MK7", "MK8", "MK9", "MK10", "MK11", "MK12", "MK13", "combined"))

pdf(paste0(Sys.Date(), "_Menopause_microbiome_FMT_MK_quantification_boxplot_stool.pdf"), width = 10, height = 5)
ggplot(vk2.quant.stool.long, aes(x=MK, y=Value, fill=FMT)) + 
  geom_boxplot(outlier.shape = NA, position=position_dodge(0.8)) +
  geom_jitter(color="black", size=1, alpha=0.9, position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75)) +
  scale_fill_manual(values=c("FMT_YF"="deeppink1", "FMT_OF"="deeppink4")) +
  theme_bw() +
  xlab("VKs") + 
  ylab("Quantification (pmol/g)") +
  ggtitle("VK quantification - Stool")
dev.off()

wilcox.test(vk2.quant.stool$MK5[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK5[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 1
wilcox.test(vk2.quant.stool$MK6[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK6[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.1905
wilcox.test(vk2.quant.stool$MK7[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK7[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.2683
wilcox.test(vk2.quant.stool$MK8[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK8[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 1
wilcox.test(vk2.quant.stool$MK9[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK9[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.7302
wilcox.test(vk2.quant.stool$MK10[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK10[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.4127
wilcox.test(vk2.quant.stool$MK11[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK11[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.4127
wilcox.test(vk2.quant.stool$MK12[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK12[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.5556
wilcox.test(vk2.quant.stool$MK13[vk2.quant.stool$FMT == "FMT_YF"], vk2.quant.stool$MK13[vk2.quant.stool$FMT == "FMT_OF"])    # p-value = 0.1905

my.MKs <- list ("FMT_YF"  = vk2.quant.stool$combined[vk2.quant.stool$FMT %in% "FMT_YF"],
                "FMT_EF"  = vk2.quant.stool$combined[vk2.quant.stool$FMT %in% "FMT_OF"])

wilcox.test(vk2.quant.stool$combined[vk2.quant.stool$FMT %in% "FMT_YF"], 
            vk2.quant.stool$combined[vk2.quant.stool$FMT %in% "FMT_OF"])     # p-value = 0.4127

pdf(paste0(Sys.Date(), "_Menopause_microbiome_FMT_MK_quantification_boxplot_stool_combined.pdf"), width = 10, height = 5)
boxplot(my.MKs, 
        outline = FALSE,
        ylim = c(0, 10000),
        col = c("purple", "orange"),
        las = 1,
        ylab = "Quantification (pmol/g)",
        main = "MK quantification - Stool")
beeswarm(my.MKs, pch = 16, col = "black", add = TRUE, cex = 1)
text(1.5,10,"p-value = 0.4127")
dev.off()

###########################
# Stool - FC
###########################

# Create a vector of column names to be used
column_names <- c("MK5", "MK6", "MK7", "MK8", "MK9", "MK10", "MK11", "MK12", "MK13")

# Filter FMT_YF rows and calculate the median for each column
median_YF <- vk2.quant.stool %>%
  filter(FMT == "FMT_YF") %>%
  summarise(across(all_of(column_names), median), .groups = 'drop')

# Convert to a list
median_YF_list <- as.list(unlist(median_YF))

# Start with the original dataframe
vk2.quant.stool.FC <- vk2.quant.stool

# Loop through column names and calculate fold change
for (col in column_names){
  vk2.quant.stool.FC[[paste0("fc_", col)]] <- vk2.quant.stool[[col]] / median_YF_list[[col]]
}

# Now remove the original columns to have only fold change columns
vk2.quant.stool.FC <- vk2.quant.stool.FC %>% 
  select(starts_with("fc_"))

vk2.quant.stool.FC <- vk2.quant.stool.FC %>% select(-fc_MK5)
vk2.quant.stool.FC$Mouse.ID <- vk2.quant.stool$Mouse.ID
vk2.quant.stool.FC$FMT <- vk2.quant.stool$FMT

# grouped boxplot

vk2.quant.stool.FC.long <- vk2.quant.stool.FC %>%
  pivot_longer(cols = fc_MK6:fc_MK13, names_to = "MK", values_to = "Value")

vk2.quant.stool.FC.long$FMT <- factor(vk2.quant.stool.FC.long$FMT, levels = c("FMT_YF", "FMT_OF"))
vk2.quant.stool.FC.long$MK <- factor(vk2.quant.stool.FC.long$MK, levels = c("fc_MK6", "fc_MK7", "fc_MK8", "fc_MK9", "fc_MK10", "fc_MK11", "fc_MK12", "fc_MK13"))

pdf(paste0(Sys.Date(), "_Menopause_microbiome_FMT_VK_quantification_boxplot_stool_FC.pdf"), width = 10, height = 5)
ggplot(vk2.quant.stool.FC.long, aes(x=MK, y=Value, fill=FMT)) + 
  geom_boxplot(outlier.shape = NA, position=position_dodge(0.8)) +
  geom_jitter(color="black", size=1, alpha=0.9, position=position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75)) +
  scale_fill_manual(values=c("FMT_YF"="purple", "FMT_OF"="orange")) +
  theme_bw() +
  xlab("VKs") + 
  ylab("Quantification (FC)") +
  ggtitle("VK quantification - Stool FC")
dev.off()

