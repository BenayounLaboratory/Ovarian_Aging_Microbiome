setwd("/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Github/Ovarian_health_index/")

library('beeswarm')

###################################
# Menopause-microbiome project
# Calculate ovarian health index scores
# Parameters: AMH, FSH, INHBA
#             Follicle count (combined)
# Cohorts: AC females, VCD, FMT
###################################

####################################
# 1. Import data files
###################################

# AC - females

my.hormone.data.AC.females <- read.table("./Input_files/MM_AC_cohort_females_hormone_data.txt", header= TRUE, sep = "\t")
my.follicle.data.AC.females <- read.table("./Input_files/MM_AC_cohort_females_ovarian_follicle_count_median.txt", header= TRUE, sep = "\t")

my.follicle.data.AC.females$combined <- rowSums(my.follicle.data.AC.females[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

# VCD cohort

my.hormone.data.VCD <- read.table("./Input_files/MM_VCD_cohort_hormone_data.txt", header= TRUE, sep = "\t")
my.follicle.data.VCD <- read.table("./Input_files/MM_VCD_cohort_ovarian_follicle_count_median.txt", header= TRUE, sep = "\t")

my.follicle.data.VCD$combined <- rowSums(my.follicle.data.VCD[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

# FMT cohort

my.hormone.data.FMT <- read.table("./Input_files/MM_FMT_cohort_hormone_data.txt", header= TRUE, sep = "\t")
my.follicle.data.FMT <- read.table("./Input_files/MM_FMT_cohort_ovarian_follicle_count_median.txt", header= TRUE, sep = "\t")

my.follicle.data.FMT$combined <- rowSums(my.follicle.data.FMT[, c("Primordial", "Primary", "Secondary", "Antral", "CL")], na.rm=TRUE)

####################################
# 2. Assign scores based on quartiles
###################################

# AC cohort - females

my.data.AC.females <- my.follicle.data.AC.females[,c("Mouse_ID", "Age", "combined")]

# For AMH, higher is better
my.data.AC.females$AMH_score <- as.numeric(cut(my.hormone.data.AC.females$AMH, breaks=quantile(my.hormone.data.AC.females$AMH), labels=c(1,2,3,4), include.lowest=TRUE))

# For FSH and INHBA, lower is better
my.data.AC.females$FSH_score <- as.numeric(cut(my.hormone.data.AC.females$FSH, breaks=quantile(my.hormone.data.AC.females$FSH), labels=c(4,3,2,1), include.lowest=TRUE))
my.data.AC.females$INHBA_score <- as.numeric(cut(my.hormone.data.AC.females$INHBA, breaks=quantile(my.hormone.data.AC.females$INHBA), labels=c(4,3,2,1), include.lowest=TRUE))

# For Follicle count, higher is better
# Manually setting break points - due to too many zero values
break_points = c(-Inf, 0, max(my.follicle.data.AC.females$combined)/3, 2*max(my.follicle.data.AC.females$combined)/3, Inf)

# Use these breaks in the cut function
my.data.AC.females$follicleCounts <- as.numeric(cut(my.follicle.data.AC.females$combined,
                                                breaks = break_points,
                                                labels = c(1, 2, 3, 4),
                                                include.lowest = TRUE))

# Compute the average score
my.data.AC.females$Ovarian_Health_Index <- rowMeans(my.data.AC.females[, c("AMH_score", "FSH_score", "INHBA_score", "follicleCounts")], na.rm=TRUE)

# If you wish to scale to a 0-100 scale
my.data.AC.females$Ovarian_Health_Index_100 <- (my.data.AC.females$Ovarian_Health_Index / 4) * 100

# View the results
my.data.AC.females[c("Mouse_ID", "Age", "Ovarian_Health_Index", "Ovarian_Health_Index_100")]

write.table(my.data.AC.females, file = "./MM_Ovarian_health_score_AC_cohort_females_results.txt", row.names = FALSE, quote = FALSE)


# VCD cohort

my.data.VCD <- my.follicle.data.VCD[,c("Mouse_ID", "Treatment", "combined")]

# For AMH, higher is better
my.data.VCD$AMH_score <- as.numeric(cut(my.hormone.data.VCD$AMH, breaks=quantile(my.hormone.data.VCD$AMH), labels=c(1,2,3,4), include.lowest=TRUE))

# For FSH and INHBA, lower is better
my.data.VCD$FSH_score <- as.numeric(cut(my.hormone.data.VCD$FSH, breaks=quantile(my.hormone.data.VCD$FSH), labels=c(4,3,2,1), include.lowest=TRUE))
my.data.VCD$INHBA_score <- as.numeric(cut(my.hormone.data.VCD$INHBA, breaks=quantile(my.hormone.data.VCD$INHBA), labels=c(4,3,2,1), include.lowest=TRUE))

# For Follicle count, higher is better
my.data.VCD$follicleCounts <- as.numeric(cut(my.follicle.data.VCD$combined, breaks=quantile(my.follicle.data.VCD$combined), labels=c(1,2,3,4), include.lowest=TRUE))

# Compute the average score
my.data.VCD$Ovarian_Health_Index <- rowMeans(my.data.VCD[, c("AMH_score", "FSH_score", "INHBA_score", "follicleCounts")], na.rm=TRUE)

# If you wish to scale to a 0-100 scale
my.data.VCD$Ovarian_Health_Index_100 <- (my.data.VCD$Ovarian_Health_Index / 4) * 100

# View the results
my.data.VCD[c("Mouse_ID", "Treatment", "Ovarian_Health_Index", "Ovarian_Health_Index_100")]

write.table(my.data.VCD, file = "./MM_Ovarian_health_score_VCD_cohort_results.txt", row.names = FALSE, quote = FALSE)


# FMT cohort

my.data.FMT <- my.follicle.data.FMT[,c("MouseID", "FMT", "combined")]

# For AMH, higher is better
my.data.FMT$AMH_score <- as.numeric(cut(my.hormone.data.FMT$AMH, breaks=quantile(my.hormone.data.FMT$AMH), labels=c(1,2,3,4), include.lowest=TRUE))

# For FSH and INHBA, lower is better
my.data.FMT$FSH_score <- as.numeric(cut(my.hormone.data.FMT$FSH, breaks=quantile(my.hormone.data.FMT$FSH), labels=c(4,3,2,1), include.lowest=TRUE))
my.data.FMT$INHBA_score <- as.numeric(cut(my.hormone.data.FMT$INHBA, breaks=quantile(my.hormone.data.FMT$INHBA), labels=c(4,3,2,1), include.lowest=TRUE))

# For Follicle count, higher is better
my.data.FMT$follicleCounts <- as.numeric(cut(my.follicle.data.FMT$combined, breaks=quantile(my.follicle.data.FMT$combined), labels=c(1,2,3,4), include.lowest=TRUE))

# Compute the average score
my.data.FMT$Ovarian_Health_Index <- rowMeans(my.data.FMT[, c("AMH_score", "FSH_score", "INHBA_score", "follicleCounts")], na.rm=TRUE)

# If you wish to scale to a 0-100 scale
my.data.FMT$Ovarian_Health_Index_100 <- (my.data.FMT$Ovarian_Health_Index / 4) * 100

# View the results
my.data.FMT[c("MouseID", "FMT", "Ovarian_Health_Index", "Ovarian_Health_Index_100")]

write.table(my.data.FMT, file = "./MM_Ovarian_health_score_FMT_cohort_results.txt", row.names = FALSE, quote = FALSE)

####################################
# 3. Plot boxplots
###################################

# AC cohort - females

my.index.AC.females <- list ("YF"  = my.data.AC.females$Ovarian_Health_Index_100[my.data.AC.females$Age %in% "YF"],
                             "EF"  = my.data.AC.females$Ovarian_Health_Index_100[my.data.AC.females$Age %in% "EF"])

wilcox.test(my.data.AC.females$Ovarian_Health_Index_100[my.data.AC.females$Age == "YF"], 
            my.data.AC.females$Ovarian_Health_Index_100[my.data.AC.females$Age == "EF"])     # p-value = 0.01091

pdf(paste(Sys.Date(),"MM_AC_cohort_females_ovarian_index_score_boxplot.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.index.AC.females, 
        outline = FALSE,
        ylim = c(0, 100),
        col = c("deeppink1", "deeppink4"),
        las = 1,
        ylab = "Ovarian health index",
        main = "AC females - Ovarian health index")
beeswarm(my.index.AC.females, pch = 16, col = "black", add = TRUE, cex = 1)
text(1.5,100,"p-value = 0.01091")
dev.off()


# VCD cohort

my.index.VCD <- list ("CTL"  = my.data.VCD$Ovarian_Health_Index_100[my.data.VCD$Treatment %in% "CTL"],
                      "VCD"  = my.data.VCD$Ovarian_Health_Index_100[my.data.VCD$Treatment %in% "VCD"])

wilcox.test(my.data.VCD$Ovarian_Health_Index_100[my.data.VCD$Treatment == "CTL"], my.data.VCD$Ovarian_Health_Index_100[my.data.VCD$Treatment == "VCD"])     # p-value = 0.01042

pdf(paste(Sys.Date(),"MM_VCD_cohort_ovarian_index_score_boxplot.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.index.VCD, 
        outline = FALSE,
        ylim = c(0, 100),
        col = c("green", "yellow"),
        las = 1,
        ylab = "Ovarian health index",
        main = "VCD cohort - Ovarian health index")
beeswarm(my.index.VCD, pch = 16, col = "black", add = TRUE, cex = 1)
text(1.5,100,"p-value = 0.01042")
dev.off()


# FMT cohort

my.index.FMT <- list ("FMT_YF"  = my.data.FMT$Ovarian_Health_Index_100[my.data.FMT$FMT %in% "YF_FMT"],
                      "FMT_EF"  = my.data.FMT$Ovarian_Health_Index_100[my.data.FMT$FMT %in% "EF_FMT"])

wilcox.test(my.data.FMT$Ovarian_Health_Index_100[my.data.FMT$FMT == "YF_FMT"], my.data.FMT$Ovarian_Health_Index_100[my.data.FMT$FMT == "EF_FMT"])     # p-value = 0.04449

pdf(paste(Sys.Date(),"MM_FMT_cohort_ovarian_index_score_boxplot.pdf", sep = "_"), width = 5, height = 7)
boxplot(my.index.FMT, 
        outline = FALSE,
        ylim = c(0, 100),
        col = c("purple", "orange"),
        las = 1,
        ylab = "Ovarian health index",
        main = "FMT cohort - Ovarian health index")
beeswarm(my.index.FMT, pch = 16, col = "black", add = TRUE, cex = 1)
text(1.5,100,"p-value = 0.04449")
dev.off()
