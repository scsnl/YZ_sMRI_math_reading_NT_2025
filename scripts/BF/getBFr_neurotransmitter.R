# -------------------------------------------------------------------------
# Replication Bayes Factor (BF10) Analysis
# Author: Yuan Zhang
# Date: 2025-07-25 
#
# This script computes replication Bayes Factors between CMI and Stanford
# datasets for:
#   - Math-related GMV maps
#   - Reading-related GMV maps
#
# For each receptor, BF10 values are computed for both directions:
#   - Using CMI as original and Stanford as replication (this is our focus)
#   - Using Stanford as original and CMI as replication (out of curiosity)
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(BayesFactor)

# -------------------------------------------------------------------------
# Load neurotransmitter receptor data
# -------------------------------------------------------------------------
fname_receptors <- 'data/neurotransmitter/receptor_data_bn246.csv'
receptors <- read.csv(fname_receptors, header = FALSE)

# Define receptor names (columns of receptor_data_bn246.csv)
receptor_names <- c(
  '5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
  'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
  'NMDA', 'VAChT'
)

# -------------------------------------------------------------------------
# Function to compute replication Bayes Factors for one modality (math/reading)
# -------------------------------------------------------------------------
compute_replication_bf <- function(fname_cmi, fname_su, bmap_name, output_prefix) {
  # Load brain-behavior maps for CMI and Stanford
  bmaps_cmi <- read.csv(fname_cmi)
  bmaps_su <- read.csv(fname_su)
  
  # Combine GMV maps with receptor data (first 218 rows only)
  df_cmi <- data.frame(bmaps_cmi, receptors[1:218, ])
  df_su <- data.frame(bmaps_su, receptors[1:218, ])
  
  # Assign receptor column names
  colnames(df_cmi)[5:ncol(df_cmi)] <- receptor_names
  colnames(df_cmi) <- make.names(colnames(df_cmi))
  colnames(df_su) <- colnames(df_cmi)  # Match column names for replication dataset
  
  # Age bins (columns 1–4 in the GMV weight CSV)
  agebins <- c("Full", "7-8.9y", "9-10.9y", "11-13.9y")
  
  # Iterate through each age bin
  for (k in 1:length(agebins)) {
    cat("Processing", bmap_name, "-", agebins[k], "\n")
    
    # Output file for this age bin
    outputf <- paste0(output_prefix, 'replication_bf10_', bmap_name, '_receptors_', agebins[k], '.csv')
    
    # Initialize BF10 ratios
    bf10_ratio1 <- c()  # CMI as original, SU as replication
    bf10_ratio2 <- c()  # SU as original, CMI as replication
    
    # Compute BF10 for each receptor
    for (i in 5:ncol(df_cmi)) {
      formula <- as.formula(paste(colnames(df_cmi)[k], " ~ ", colnames(df_cmi)[i], sep = ""))
      
      # BF for original dataset (CMI)
      BF1 <- lmBF(formula, data = df_cmi)
      
      # BF for replication dataset (Stanford)
      BF2 <- lmBF(formula, data = df_su)
      
      # BF for combined dataset
      BF_total <- lmBF(formula, data = rbind(df_cmi, df_su))
      
      # Compute ratios
      ratio1 <- extractBF(BF_total, logbf = FALSE, onlybf = TRUE) / extractBF(BF1, logbf = FALSE, onlybf = TRUE)
      ratio2 <- extractBF(BF_total, logbf = FALSE, onlybf = TRUE) / extractBF(BF2, logbf = FALSE, onlybf = TRUE)
      
      bf10_ratio1 <- c(bf10_ratio1, ratio1)
      bf10_ratio2 <- c(bf10_ratio2, ratio2)
    }
    
    # Create results dataframe
    res <- data.frame(receptor = receptor_names, BF10_cmi_su = bf10_ratio1, BF10_su_cmi = bf10_ratio2)
    res$bmap <- bmap_name
    res$agbin <- agebins[k]
    
    # Reorder rows by a custom receptor order
    idx <- c(9:11, 15, 18, 12, 1:7, 14, 19, 17, 8, 13, 16)
    res <- res[idx, c(4, 1, 2, 3, 5)]
    
    # Save results
    write.csv(res, file = outputf, row.names = FALSE)
    cat("Replication BF results saved to:", outputf, "\n")
  }
}

# -------------------------------------------------------------------------
# Run replication BF analysis for Math
# -------------------------------------------------------------------------
compute_replication_bf(
  fname_cmi = 'results/age_analysis/cmi/gmv_weights_agebin_math.csv',
  fname_su = 'results/age_analysis/stanford/gmv_weights_agebin_math.csv',
  bmap_name = 'math_gmv',
  output_prefix = 'results/neurotransmitter/stanford/'
)

# -------------------------------------------------------------------------
# Run replication BF analysis for Reading
# -------------------------------------------------------------------------
compute_replication_bf(
  fname_cmi = 'results/age_analysis/cmi/gmv_weights_agebin_read.csv',
  fname_su = 'results/age_analysis/stanford/gmv_weights_agebin_read.csv',
  bmap_name = 'read_gmv',
  output_prefix = 'results/neurotransmitter/stanford/'
)
