setwd("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Serum_hormone_quantification_FMT")

library(tidyverse)
library(ggplot2)

##############################################
# 1. Load data and assess relationship
##############################################

# Import data
my.data <- read.table("/Users/minhookim/Dropbox/Benayoun_lab/Menopause_microbiome_project/For_revision/Data_analysis/Serum_hormone_quantification_FMT/FSH_LH_Multiplex_Utrasensitive.txt", header = TRUE, sep = "\t")

# Scatterplot for FSH
ggplot(my.data, aes(x = FSH_Multiplex, y = FSH_US)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "FSH: Multiplex vs. US", x = "Multiplex FSH", y = "US FSH")

##############################################
# 2. Fit models for normalization
##############################################

# Linear model
fsh_lm <- lm(FSH_Multiplex ~ FSH_US, data = my.data)

# Polynomial model 
fsh_poly <- lm(FSH_Multiplex ~ poly(FSH_US, 2), data = my.data)

# Log transformation Model
fsh_log <- lm(FSH_Multiplex ~ log(FSH_US + 1), data = my.data)

##############################################
# 3. Evaluate models
##############################################

# Compare models using R², adjusted R², and AIC
model_metrics <- function(model) {
  data.frame(
    R_squared = summary(model)$r.squared,
    Adjusted_R_squared = summary(model)$adj.r.squared,
    AIC = AIC(model)
  )
}

# FSH Models
fsh_models <- list(Linear = fsh_lm, Polynomial = fsh_poly, Log = fsh_log)
fsh_metrics <- bind_rows(lapply(fsh_models, model_metrics), .id = "Model")

print(fsh_metrics)
#        Model R_squared Adjusted_R_squared      AIC
# 1     Linear 0.7278230          0.7187504 304.8829
# 2 Polynomial 0.9001780          0.8932937 274.7849
# 3        Log 0.4933543          0.4764661 324.7664

##############################################
# 4. Use polynomial model 
##############################################

best_fsh_model <- fsh_poly 

# Predict Multiplex-equivalent values for future US data
my.data$FSH_Normalized <- predict(best_fsh_model, newdata = my.data)

ggplot(my.data, aes(x = FSH_Multiplex, y = FSH_Normalized)) +
  geom_point(color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "FSH Multiplex vs. FSH Normalized",
    x = "FSH Multiplex",
    y = "FSH Normalized"
  ) +
  theme_minimal()

# Calculate Percentage Difference and Error
my.data <- my.data %>%
  mutate(
    FSH_Percent_Diff = abs(FSH_Multiplex - FSH_Normalized) / FSH_Multiplex * 100,
    FSH_Absolute_Error = abs(FSH_Multiplex - FSH_Normalized)
  )

# Summarize errors
error_summary <- my.data %>%
  summarize(
    Mean_Percent_Diff = mean(FSH_Percent_Diff, na.rm = TRUE),
    Mean_Absolute_Error = mean(FSH_Absolute_Error, na.rm = TRUE),
    SD_Percent_Diff = sd(FSH_Percent_Diff, na.rm = TRUE),
    SD_Absolute_Error = sd(FSH_Absolute_Error, na.rm = TRUE)
  )

print(error_summary)
#   Mean_Percent_Diff Mean_Absolute_Error SD_Percent_Diff SD_Absolute_Error
# 1          414.0911            9.336833        1024.237          12.74292

ggplot(my.data, aes(x = Group, y = FSH_Percent_Diff, color = Group)) +
  geom_point(size = 3) +
  labs(
    title = "Percentage Difference Between FSH Multiplex and Normalized",
    x = "Sample ID",
    y = "Percent Difference (%)"
  ) +
  theme_minimal()


# Save model

save(fsh_poly, file = "UVA_FSH_Multiplex_US_normalization_model.RData")

################################################################################
# Save session info
sink(file = paste0(Sys.Date(), "_MM_UVA_FSH_multiplex_vs_US_normalization_Session_Info.txt"))
sessionInfo()
sink()
