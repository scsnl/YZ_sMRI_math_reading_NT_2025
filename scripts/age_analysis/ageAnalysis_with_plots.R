##########################################################################
# Age Analysis of CCA GMV Weights for CMI-HBN and Stanford Cohorts
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script performs the following:
#   1. Loads math and reading CCA results for both CMI-HBN and Stanford.
#   2. Extracts the brain-side CCA scores (U) for math and reading.
#   3. Examines correlations between CCA brain scores and age.
#   4. Computes GMV weight maps (full sample and age bins).
#   5. Saves GMV weight maps to CSV files for both math and reading.
#   6. Creates correlation plots comparing GMV weight maps across age bins.
###########################################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ggplot2)
library(Hmisc)
library(corrplot)

# ------------------------------------------------------------------
# 2. Define Datasets
# ------------------------------------------------------------------
datasets <- list(
  CMI = list(
    math_path = "results/cca/cmi/wholebrain_cca_cmi_math/",
    math_file = "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    read_path = "results/cca/cmi/wholebrain_cca_cmi_reading/",
    read_file = "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_dir = "results/age_analysis/cmi/"
  ),
  Stanford = list(
    math_path = "results/cca/stanford/wholebrain_cca_stanford_math/",
    math_file = "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    read_path = "results/cca/stanford/wholebrain_cca_stanford_reading/",
    read_file = "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_dir = "results/age_analysis/stanford/"
  )
)

# Ensure output directories exist
for (ds in names(datasets)) {
  dir.create(datasets[[ds]]$output_dir, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------------
# 3. Function to Run Age Analysis
# ------------------------------------------------------------------
analyze_age <- function(math_path, math_file, read_path, read_file, output_dir, dataset_name) {
  cat("\n===================================\n")
  cat("Processing Dataset:", dataset_name, "\n")
  cat("===================================\n")
  
  # --- Load Math CCA Results ---
  load(paste0(math_path, math_file))
  Umath <- res$Cx[, 2]   # Brain-side CCA mode 2 scores
  gmv <- data$brain      # GMV data (nsub x nroi)
  beh <- data$beh        # Behavioral data
  
  # --- Load Reading CCA Results ---
  load(paste0(read_path, read_file))
  Uread <- res$Cx[, 2]   # Brain-side CCA mode 2 scores
  
  # --- Correlation with Age ---
  cat("\nCorrelation with Age (Math):\n")
  print(cor.test(-1 * Umath, beh$age))
  
  cat("\nCorrelation with Age (Reading):\n")
  print(cor.test(-1 * Uread, beh$age))
  
  # ------------------------------------------------------------------
  # 4. Compute GMV Weight Maps (Full Sample & Age Bins)
  # ------------------------------------------------------------------
  xloading_math_full <- cor(gmv, -1 * Umath)
  xloading_read_full <- cor(gmv, -1 * Uread)
  
  # Age bins: 7–8.9, 9–10.9, 11–13.9
  bin_labels <- c("7-8.9y", "9-10.9y", "11-13.9y")
  
  # Bin 1: 7–8.9
  idx <- which(beh$age < 9)
  xloading_math_bin1 <- cor(gmv[idx, ], -1 * Umath[idx])
  xloading_read_bin1 <- cor(gmv[idx, ], -1 * Uread[idx])
  
  # Bin 2: 9–10.9
  idx <- which(beh$age >= 9 & beh$age < 11)
  xloading_math_bin2 <- cor(gmv[idx, ], -1 * Umath[idx])
  xloading_read_bin2 <- cor(gmv[idx, ], -1 * Uread[idx])
  
  # Bin 3: 11–13.9
  idx <- which(beh$age >= 11)
  xloading_math_bin3 <- cor(gmv[idx, ], -1 * Umath[idx])
  xloading_read_bin3 <- cor(gmv[idx, ], -1 * Uread[idx])
  
  # Combine into dataframes
  df_math <- data.frame(
    xloading_math_full, xloading_math_bin1, xloading_math_bin2, xloading_math_bin3
  )
  colnames(df_math) <- c("Full", bin_labels)
  
  df_read <- data.frame(
    xloading_read_full, xloading_read_bin1, xloading_read_bin2, xloading_read_bin3
  )
  colnames(df_read) <- c("Full", bin_labels)
  
  # Save to CSV
  write.csv(df_math, paste0(output_dir, "gmv_weights_agebin_math.csv"), row.names = FALSE)
  write.csv(df_read, paste0(output_dir, "gmv_weights_agebin_read.csv"), row.names = FALSE)
  
  # ------------------------------------------------------------------
  # 5. Correlation Plots for Math & Reading GMV Weights
  # ------------------------------------------------------------------
  plot_corr <- function(df, output_file) {
    stats <- rcorr(as.matrix(df))  # Pearson correlation matrix
    cor_matrix <- stats$r
    colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(df)
    
    # Save source data for the plotted correlation matrix
    source_data_file <- sub("\\.eps$", "_source_data_correlation_matrix.csv", output_file)
    write.csv(cor_matrix, file = source_data_file, row.names = TRUE)
    
    # also save p-value matrix, if useful for checking
    pval_file <- sub("\\.eps$", "_source_data_pvalue_matrix.csv", output_file)
    write.csv(stats$P, file = pval_file, row.names = TRUE)
    
    # Custom color palette
    my_col <- colorRampPalette(c("orange", "tomato"))(200)
    
    corrplot(
      cor_matrix, method = "color", col = my_col, type = "upper", diag = FALSE,
      addCoef.col = "white", number.cex = 0.85,
      tl.col = "black", tl.srt = 45, cl.pos = "n",
      is.corr = FALSE, cl.lim = c(0, 1), mar = c(0, 0, 1, 0)
    )
    dev.copy2eps(file = output_file, width = 5, height = 5)
  }
  
  plot_corr(df_math, paste0(output_dir, "gmv_weights_correlation_math.eps"))
  plot_corr(df_read, paste0(output_dir, "gmv_weights_correlation_read.eps"))
}

# ------------------------------------------------------------------
# 6. Run Analysis for Both Datasets
# ------------------------------------------------------------------
for (ds in names(datasets)) {
  analyze_age(
    math_path = datasets[[ds]]$math_path,
    math_file = datasets[[ds]]$math_file,
    read_path = datasets[[ds]]$read_path,
    read_file = datasets[[ds]]$read_file,
    output_dir = datasets[[ds]]$output_dir,
    dataset_name = ds
  )
}
