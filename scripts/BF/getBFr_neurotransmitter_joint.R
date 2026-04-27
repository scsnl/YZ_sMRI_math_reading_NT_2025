# -------------------------------------------------------------------------
# Replication Bayes Factor (BF10) Analysis for current domain-specific results
# Author: Yuan Zhang
# Date: 2026-04-13
#
# This script computes replication Bayes Factors between CMI and Stanford
# datasets for:
#   - joint CCA Mode 2 GMV map
#
# Brain map source:
#   *_coef.csv
# Weight column used:
#   xloading_m2
# Sign is flipped by multiplying by -1 before analysis
#
# For each receptor, BF10 values are computed for both directions:
#   - Using CMI as original and Stanford as replication (our focus)
#   - Using Stanford as original and CMI as replication
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(BayesFactor)

# -------------------------------------------------------------------------
# Load neurotransmitter receptor data
# -------------------------------------------------------------------------
fname_receptors <- "data/neurotransmitter/receptor_data_bn246.csv"
receptors <- read.csv(fname_receptors, header = FALSE)

# Define receptor names
receptor_names <- c(
  "5HT1a", "5HT1b", "5HT2a", "5HT4", "5HT6", "5HTT", "A4B2", "CB1",
  "D1", "D2", "DAT", "GABAa", "H3", "M1", "mGluR5", "MOR", "NET",
  "NMDA", "VAChT"
)

# -------------------------------------------------------------------------
# Function to compute replication Bayes Factors for one domain-specific map
# -------------------------------------------------------------------------
compute_replication_bf_mode2 <- function(fname_cmi, fname_su, bmap_name, output_file,
                                         weight_col = "xloading_m2",
                                         flip_sign = TRUE) {
  # Load coefficient CSVs
  coef_cmi <- read.csv(fname_cmi)
  coef_su  <- read.csv(fname_su)
  
  if (!(weight_col %in% colnames(coef_cmi))) {
    stop(paste(weight_col, "not found in", fname_cmi))
  }
  if (!(weight_col %in% colnames(coef_su))) {
    stop(paste(weight_col, "not found in", fname_su))
  }
  
  # Extract mode-2 brain maps
  bmap_cmi <- coef_cmi[[weight_col]][1:218]
  bmap_su  <- coef_su[[weight_col]][1:218]
  
  # Flip sign if requested
  if (flip_sign) {
    bmap_cmi <- -1 * bmap_cmi
    bmap_su  <- -1 * bmap_su
  }
  
  # Build dataframes
  df_cmi <- data.frame(
    bmap = bmap_cmi,
    receptors[1:218, ]
  )
  
  df_su <- data.frame(
    bmap = bmap_su,
    receptors[1:218, ]
  )
  
  # Assign receptor names
  colnames(df_cmi)[2:ncol(df_cmi)] <- receptor_names
  colnames(df_su)[2:ncol(df_su)] <- receptor_names
  
  # Ensure valid column names
  colnames(df_cmi) <- make.names(colnames(df_cmi))
  colnames(df_su)  <- make.names(colnames(df_su))
  
  # Compute BF10 for each receptor
  bf10_ratio1 <- c()  # CMI as original, Stanford as replication
  bf10_ratio2 <- c()  # Stanford as original, CMI as replication
  
  for (i in 2:ncol(df_cmi)) {
    formula <- as.formula(paste("bmap ~", colnames(df_cmi)[i]))
    
    # BF for original dataset (CMI)
    BF1 <- lmBF(formula, data = df_cmi)
    
    # BF for replication dataset (Stanford)
    BF2 <- lmBF(formula, data = df_su)
    
    # BF for combined dataset
    BF_total <- lmBF(formula, data = rbind(df_cmi, df_su))
    
    # Compute replication BF ratios
    ratio1 <- extractBF(BF_total, logbf = FALSE, onlybf = TRUE) /
      extractBF(BF1, logbf = FALSE, onlybf = TRUE)
    
    ratio2 <- extractBF(BF_total, logbf = FALSE, onlybf = TRUE) /
      extractBF(BF2, logbf = FALSE, onlybf = TRUE)
    
    bf10_ratio1 <- c(bf10_ratio1, ratio1)
    bf10_ratio2 <- c(bf10_ratio2, ratio2)
  }
  
  # Create results dataframe
  res <- data.frame(
    bmap = bmap_name,
    receptor = receptor_names,
    BF10_cmi_su = bf10_ratio1,
    BF10_su_cmi = bf10_ratio2
  )
  
  # Reorder rows to match previous output style
  idx <- c(9:11, 15, 18, 12, 1:7, 14, 19, 17, 8, 13, 16)
  res <- res[idx, c("bmap", "receptor", "BF10_cmi_su", "BF10_su_cmi")]
  
  # Save results
  write.csv(res, file = output_file, row.names = FALSE)
  cat("Replication BF results saved to:", output_file, "\n")
}

# -------------------------------------------------------------------------
# Run replication BF analysis for shared Mode 2
# -------------------------------------------------------------------------
compute_replication_bf_mode2(
  fname_cmi = "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_coef.csv",
  fname_su  = "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined/CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_coef.csv",
  bmap_name = "shared_mode2",
  output_file = "results/neurotransmitter/stanford/shared/replication_bf10_shared_mode2_receptors.csv"
)
