options(stringsAsFactors = F)

library(ggplot2)
library(dplyr)
library('beeswarm')
library('pheatmap')
library('grDevices')
library(bitops)
library("forcats")
library(survival)  # survival_3.3-1

###################################
# Menopause-microbiome project
# Plot fertility data - FMT cohort
###################################

###################################
# Import files
###################################

my.fertility.assay <- read.table("./FMT_cohort_fertility_data_combined.txt", header = TRUE)

my.fertility.assay$FMT <- factor(my.fertility.assay$FMT , levels=c("FMT_YF", "FMT_EF"))

###################################
# Boxplots for fertility assay
###################################

my.fertility <- list ("FMT_YF"  = my.fertility.assay$Num_pups[my.fertility.assay$FMT %in% "FMT_YF"],
                      "FMT_EF"  = my.fertility.assay$Num_pups[my.fertility.assay$FMT %in% "FMT_EF"])

wilcox.test(my.fertility.assay$Num_pups[my.fertility.assay$FMT %in% "FMT_YF"], my.fertility.assay$Num_pups[my.fertility.assay$FMT %in% "FMT_EF"]) #p-value = 0.1726

# Number of pups

pdf(paste(Sys.Date(),"MM_FMT_cohort_fertility_assay_pup_count.pdf", sep = "_"), width = 6, height = 7)
boxplot(Num_pups ~ FMT, data = my.fertility.assay,
        main = "Fertility assay combined",
        ylim = c(0, 11),
        col = c("purple", "orange"),
        xlab = "FMT",
        ylab = "Litter count")
beeswarm(Num_pups ~ FMT, data = my.fertility.assay, pch = 16, add = TRUE)
text(1.5,11,"p-value = 0.1726")
dev.off()

# Latency

fit.latency <- survfit(Surv(Latency, Status) ~ FMT, data=my.fertility.assay)

#running lg rank test
surv.test <- survdiff(Surv(Latency, Status) ~ FMT, data=my.fertility.assay)
surv.test
# survdiff(formula = Surv(Latency, Status) ~ FMT, data = my.fertility.assay)

# N Observed Expected (O-E)^2/E (O-E)^2/V
# FMT=FMT_YF 18       18     24.3      1.63      9.64
# FMT=FMT_EF 18       18     11.7      3.38      9.64

# Chisq= 9.6  on 1 degrees of freedom, p= 0.002 

pdf(paste0(Sys.Date(),"_MM_FMT_cohort_reproductive_success_span.pdf"), height = 4.5, width = 5)
plot(fit.latency, col= c("purple","orange"),  lwd=2, mark.time=TRUE,
     xlab="Days since male pairing", ylab="% nulliparous", main = "Reproductive success", xlim = c(15,30), ylim = c(0, 1.1))
legend("topright", 
       c("FMT-YF (n = 16)", "FMT-EF (n = 16)"),
       col= c("purple","orange"), lwd=2, bty='n')
text(22.5, 1.1, "p = 0.002")
dev.off()

#######################
sink(file = paste(Sys.Date(),"MM_MM_cohort_fertility_data_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()
