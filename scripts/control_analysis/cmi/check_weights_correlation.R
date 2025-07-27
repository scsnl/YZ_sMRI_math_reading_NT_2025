############################################################
# Compare CCA Mode Scores Across Different Models
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script loads canonical variate scores (U, V) from 
# multiple CCA models (math and reading) for the CMI cohort, 
# including:
#   - Original models 
#   - IQ-controlled models
#   - SES-controlled models
#
# It then compares the mode scores between the original, 
# IQ-controlled, and SES-controlled models using Pearson 
# correlation tests.
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

full_id <- data$beh$oak_id      # Participant IDs for the full dataset
Umath   <- res$Cx[, 2]          # Second canonical variate (brain-side) for math
Vmath   <- res$Cy[, 2]          # Second canonical variate (behavior-side) for math

#################################
## Reading - Original CCA Model ##
#################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread <- res$Cx[, 2]            # Second canonical variate (brain-side) for reading
Vread <- res$Cy[, 2]            # Second canonical variate (behavior-side) for reading

#####################################
## Math - IQ-Controlled CCA Model  ##
#####################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math_ctrIQ/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Umath_ctrIQ <- res$Cx[, 2]      # Math brain-side score (IQ-controlled)
Vmath_ctrIQ <- res$Cy[, 2]      # Math behavior-side score (IQ-controlled)

#######################################
## Reading - IQ-Controlled CCA Model ##
#######################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading_ctrIQ/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread_ctrIQ <- res$Cx[, 2]      # Reading brain-side score (IQ-controlled)
Vread_ctrIQ <- res$Cy[, 2]      # Reading behavior-side score (IQ-controlled)

#####################################
## Math - SES-Controlled CCA Model ##
#####################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math_ctrSES/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

ses_id <- data$beh$oak_id       # Participant IDs for SES dataset
Umath_ctrSES <- res$Cx[, 2]     # Math brain-side score (SES-controlled)
Vmath_ctrSES <- res$Cy[, 2]     # Math behavior-side score (SES-controlled)

#######################################
## Reading - SES-Controlled CCA Model ##
#######################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading_ctrSES/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread_ctrSES <- res$Cx[, 2]     # Reading brain-side score (SES-controlled)
Vread_ctrSES <- res$Cy[, 2]     # Reading behavior-side score (SES-controlled)

# ------------------------------------------------------------------
# 3. Correlation Analysis
# ------------------------------------------------------------------

###############################################################
## Correlate Math/Reading Mode Scores Between Model Variants ##
###############################################################

# -------- IQ-Controlled vs Original --------
cat("\n### Correlation between Original and IQ-Controlled Models ###\n")
cor.test(Umath, Umath_ctrIQ)    # Brain-side (math)
cor.test(Vmath, Vmath_ctrIQ)    # Behavior-side (math)
cor.test(Uread, Uread_ctrIQ)    # Brain-side (reading)
cor.test(Vread, Vread_ctrIQ)    # Behavior-side (reading)

# -------- SES-Controlled vs Original --------
cat("\n### Correlation between Original and SES-Controlled Models ###\n")
idx <- which(full_id %in% ses_id)  # Match IDs (since SES model might exclude NAs)
cor.test(Umath[idx], Umath_ctrSES) # Brain-side (math)
cor.test(Vmath[idx], Vmath_ctrSES) # Behavior-side (math)
cor.test(Uread[idx], Uread_ctrSES) # Brain-side (reading)
cor.test(Vread[idx], Vread_ctrSES) # Behavior-side (reading)
