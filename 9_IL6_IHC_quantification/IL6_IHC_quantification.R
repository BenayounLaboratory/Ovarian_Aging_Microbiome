library(beeswarm)
library(outliers)

###################################
# Menopause-microbiome project - revision #2
# Generate plots for IHC data
###################################

###################################
# 1. Import quantification data
###################################

IL6_data <- read.table("./IL6_quantification.txt", header = TRUE, sep = "\t")

IL6_data$Group <- factor(IL6_data$Group, levels = c("FMT_YF", "FMT_EF"))

###################################
# 2. Calculate average of triplicates
###################################

# Calculate average percentage per Mouse_ID
IL6_ave <- aggregate(Percentage ~ Mouse_ID + Group, data = IL6_data, FUN = mean)
IL6_ave$Cohort <- ifelse(grepl("MM", IL6_ave$Mouse_ID), "Cohort1", "Cohort2")

# Wilcoxon test
wilcox_result_IL6 <- wilcox.test(IL6_ave$Percentage[IL6_ave$Group == "FMT_YF"], 
                                 IL6_ave$Percentage[IL6_ave$Group == "FMT_EF"])

# Generate boxplot
pdf(paste(Sys.Date(), "MM_FMT_cohort_IHC_IL6_quantification.pdf", sep = "_"), width = 6, height = 7)
boxplot(Percentage ~ Group, data = IL6_ave, outline = FALSE,
        main = "IHC quantification: a-IL6",
        ylim = c(0, max(IL6_ave$Percentage) * 1.1),
        col = c("purple", "orange"),
        xlab = "FMT",
        ylab = "IL6 % H DAB signal / area")
beeswarm(Percentage ~ Group, data = IL6_ave, pch = 16, add = TRUE)
text(1.5, max(IL6_ave$Percentage) * 1.05, 
     paste("p-value =", round(wilcox_result_IL6$p.value, 4)))
dev.off()

grubbs.test(IL6_ave$Percentage[IL6_ave$Group == "FMT_EF"])     # p-value = 0.001684

# Remove outlier

IL6_ave.clean <- IL6_ave[IL6_ave$Mouse_ID != "MM070", ]

# Wilcoxon test
wilcox_result_IL6 <- wilcox.test(IL6_ave.clean$Percentage[IL6_ave.clean$Group == "FMT_YF"], 
                                 IL6_ave.clean$Percentage[IL6_ave.clean$Group == "FMT_EF"])

my.pwpch <- ifelse(grepl("1", IL6_ave.clean$Cohort), 16,17)

pdf(paste(Sys.Date(), "MM_FMT_cohort_IHC_IL6_quantification_without_outlier.pdf", sep = "_"), width = 6, height = 7)
boxplot(Percentage ~ Group, data = IL6_ave.clean, outline = FALSE,
        main = "IHC quantification: a-IL6",
        ylim = c(0, 20),
        col = c("purple", "orange"),
        xlab = "FMT",
        ylab = "IL6 % H DAB signal / area")
beeswarm(Percentage ~ Group, data = IL6_ave.clean, pwpch = my.pwpch, add = TRUE)
text(1.5, 20, 
     paste("p-value =", round(wilcox_result_IL6$p.value, 4)))
dev.off()
