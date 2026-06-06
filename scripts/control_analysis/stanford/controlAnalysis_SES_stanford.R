############################################################
# Stanford SES-Control Analysis
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script:
#  1) Loads canonical variate scores (U) from math and reading 
#     CCA models for the Stanford cohort.
#  2) Computes a combined SES score from parental income and 
#     education variables.
#  3) Examines correlations between the CCA brain scores and 
#     the SES measure.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory to project root
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear the R environment
rm(list = ls())

# Load required library
library(ggplot2)

# ------------------------------------------------------------------
# 2. Load Canonical Variate Scores (CCA U modes)
# ------------------------------------------------------------------

###############################
## Math - Original CCA Model ##
###############################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

full_id <- data$beh$pid   # Participant IDs for the Stanford cohort
Umath   <- res$Cx[, 2]    # Brain-side canonical variate for math

#################################
## Reading - Original CCA Model ##
#################################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread <- res$Cx[, 2]      # Brain-side canonical variate for reading

# ------------------------------------------------------------------
# 3. Generate SES Scores
# ------------------------------------------------------------------

########################
## Generate SES Score ##
########################
beh <- data$beh   # Behavioral data frame

# Identify participants with complete income and education data for each parent
idx_parent1 <- which(!is.na(beh$income_parent1) & !is.na(beh$max_edu_parent1))
idx_parent2 <- which(!is.na(beh$income_parent2) & !is.na(beh$max_edu_parent2))
length(unique(union(idx_parent1, idx_parent2)))  # Count participants with SES info

# ------------------------------
# Z-score parental income/education
# ------------------------------
beh$z_income_parent1  <- as.vector(scale(beh$income_parent1))
beh$z_income_parent2  <- as.vector(scale(beh$income_parent2))
beh$z_max_edu_parent1 <- as.vector(scale(beh$max_edu_parent1))
beh$z_max_edu_parent2 <- as.vector(scale(beh$max_edu_parent2))

# ------------------------------
# Compute SES for each parent
# ------------------------------
beh$ses_parent1 <- (beh$z_income_parent1 + beh$z_max_edu_parent1) / 2
beh$ses_parent2 <- (beh$z_income_parent2 + beh$z_max_edu_parent2) / 2

# ------------------------------
# Combine SES across parents
# ------------------------------
beh$ses_combined <- beh$ses_parent1
for (i in 1:nrow(beh)) {
  if (!is.na(beh$ses_parent1[i]) & !is.na(beh$ses_parent2[i])) {
    # If both parents have SES info, average them
    beh$ses_combined[i] <- (beh$ses_parent1[i] + beh$ses_parent2[i]) / 2
  } else if (!is.na(beh$ses_parent2[i])) {
    # If only parent2 has SES info, use it
    beh$ses_combined[i] <- beh$ses_parent2[i]
  }
}

# Plot histogram of combined SES scores
hist(beh$ses_combined)

# ------------------------------------------------------------------
# 4. Correlation Between CCA Scores and SES
# ------------------------------------------------------------------

idx_ses <- which(!is.na(beh$ses_combined))  # Participants with valid SES

print_cor_stats <- function(x, y, label) {
  idx <- which(!is.na(x) & !is.na(y))
  ct <- cor.test(x[idx], y[idx], method = "pearson")
  
  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  cat("N =", length(idx), "\n")
  cat("Pearson r =", round(unname(ct$estimate), 4), "\n")
  cat("t(", unname(ct$parameter), ") = ", round(unname(ct$statistic), 4), "\n", sep = "")
  cat("p =", signif(ct$p.value, 4), "\n")
  
  cat("95% CI = [", round(ct$conf.int[1], 4), ", ",round(ct$conf.int[2], 4), "]\n", sep = "")
}

# Correlation between math CCA scores and SES
print_cor_stats(x = Umath, y = beh$ses_combined,
  label = "Correlation between Math CCA brain score and SES"
)

# Correlation between reading CCA scores and SES
print_cor_stats(x = Uread, y = beh$ses_combined,
  label = "Correlation between Reading CCA brain score and SES"
)
