############################################################
# Compare CCA Mode Scores Across Different Models
# Author: Yuan Zhang
# Date: 2026-03-24
#
# Description:
# This script loads canonical variate scores (U, V) from
# multiple CCA models (math and reading) for the CMI cohort,
# including:
#   - Original models
#   - IQ-controlled models
#   - SES-controlled models
#   - Site-controlled models
#
# It then compares the mode scores between the original
# and each controlled model using Pearson correlation tests.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory to project root
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear workspace
rm(list = ls())

# Load required libraries
library(ggplot2)

# ------------------------------------------------------------------
# 2. Load CCA Mode Scores
# ------------------------------------------------------------------

###############################
## Math - Original CCA Model ##
###############################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

full_id <- data$beh$oak_id
Umath   <- res$Cx[, 2]
Vmath   <- res$Cy[, 2]

#################################
## Reading - Original CCA Model ##
#################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread <- res$Cx[, 2]
Vread <- res$Cy[, 2]

#####################################
## Math - IQ-Controlled CCA Model  ##
#####################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math_ctrIQ/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Umath_ctrIQ <- res$Cx[, 2]
Vmath_ctrIQ <- res$Cy[, 2]

#######################################
## Reading - IQ-Controlled CCA Model ##
#######################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading_ctrIQ/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread_ctrIQ <- res$Cx[, 2]
Vread_ctrIQ <- res$Cy[, 2]

#####################################
## Math - SES-Controlled CCA Model ##
#####################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math_ctrSES/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

ses_id <- data$beh$oak_id
Umath_ctrSES <- res$Cx[, 2]
Vmath_ctrSES <- res$Cy[, 2]

#######################################
## Reading - SES-Controlled CCA Model ##
#######################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading_ctrSES/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread_ctrSES <- res$Cx[, 2]
Vread_ctrSES <- res$Cy[, 2]

######################################
## Math - Site-Controlled CCA Model ##
######################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math_ctrSite/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

site_id <- data$beh$oak_id
Umath_ctrSite <- res$Cx[, 2]
Vmath_ctrSite <- res$Cy[, 2]

########################################
## Reading - Site-Controlled CCA Model ##
########################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading_ctrSite/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread_ctrSite <- res$Cx[, 2]
Vread_ctrSite <- res$Cy[, 2]


# ------------------------------------------------------------------
# 3. Correlation Analysis
# ------------------------------------------------------------------

###############################################################
## Correlate Math/Reading Mode Scores Between Model Variants ##
###############################################################

# -------- IQ-Controlled vs Original --------
cat("\n### Correlation between Original and IQ-Controlled Models ###\n")
print(cor.test(Umath, Umath_ctrIQ))
print(cor.test(Vmath, Vmath_ctrIQ))
print(cor.test(Uread, Uread_ctrIQ))
print(cor.test(Vread, Vread_ctrIQ))

# -------- SES-Controlled vs Original --------
cat("\n### Correlation between Original and SES-Controlled Models ###\n")
idx <- which(full_id %in% ses_id)
print(cor.test(Umath[idx], Umath_ctrSES))
print(cor.test(Vmath[idx], Vmath_ctrSES))
print(cor.test(Uread[idx], Uread_ctrSES))
print(cor.test(Vread[idx], Vread_ctrSES))

# -------- Site-Controlled vs Original --------
cat("\n### Correlation between Original and Site-Controlled Models ###\n")
idx <- which(full_id %in% site_id)
print(cor.test(Umath[idx], Umath_ctrSite))
print(cor.test(Vmath[idx], Vmath_ctrSite))
print(cor.test(Uread[idx], Uread_ctrSite))
print(cor.test(Vread[idx], Vread_ctrSite))
