############################################################
# Compare CCA Mode Scores Across Original and IQ-Controlled Models
# (Stanford Cohort)
#
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script:
#  1. Loads the canonical variate scores (U for brain, V for behavior) 
#     from both original and IQ-controlled CCA models for the Stanford cohort.
#  2. Compares the canonical variates from original vs. IQ-controlled models
#     by calculating Pearson correlations for both math and reading tasks.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory to the project root
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear the R environment to avoid conflicts with existing objects
rm(list = ls())

# ------------------------------------------------------------------
# 2. Load CCA Mode Scores
# ------------------------------------------------------------------

###############################
## Math - Original CCA Model ##
###############################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

full_id <- data$beh$pid   # Participant IDs for math data
Umath   <- res$Cx[, 2]    # Brain-side canonical variate (2nd mode) for math
Vmath   <- res$Cy[, 2]    # Behavior-side canonical variate (2nd mode) for math

#################################
## Reading - Original CCA Model ##
#################################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread <- res$Cx[, 2]      # Brain-side canonical variate (2nd mode) for reading
Vread <- res$Cy[, 2]      # Behavior-side canonical variate (2nd mode) for reading

#####################################
## Math - IQ-Controlled CCA Model  ##
#####################################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_math_ctrIQ/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

iq_id         <- data$beh$pid  # Participant IDs for IQ-controlled math data
Umath_ctrIQ   <- res$Cx[, 2]   # Brain-side canonical variate (IQ-controlled math)
Vmath_ctrIQ   <- res$Cy[, 2]   # Behavior-side canonical variate (IQ-controlled math)

#######################################
## Reading - IQ-Controlled CCA Model ##
#######################################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_reading_ctrIQ/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep = ""))

Uread_ctrIQ <- res$Cx[, 2]     # Brain-side canonical variate (IQ-controlled reading)
Vread_ctrIQ <- res$Cy[, 2]     # Behavior-side canonical variate (IQ-controlled reading)

# ------------------------------------------------------------------
# 3. Correlation Analysis
# ------------------------------------------------------------------

###############################################################
## Check Correlation of Math/Reading Mode Scores Between Models
###############################################################

# Math: Compare original vs. IQ-controlled
cor.test(Umath, Umath_ctrIQ)   # Brain-side canonical variate correlation
cor.test(Vmath, Vmath_ctrIQ)   # Behavior-side canonical variate correlation

# Reading: Compare original vs. IQ-controlled
cor.test(Uread, Uread_ctrIQ)   # Brain-side canonical variate correlation
cor.test(Vread, Vread_ctrIQ)   # Behavior-side canonical variate correlation
