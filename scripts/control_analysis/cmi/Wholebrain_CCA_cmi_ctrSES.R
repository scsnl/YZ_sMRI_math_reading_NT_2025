############################################################
# Canonical Correlation Analysis (CCA) - CMI Cohort (SES-Controlled)
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script runs CCA analyses for both math and reading 
# tasks on the CMI cohort, controlling for age and regressing 
# out SES (Barratt_Total) from math and reading measures. 
# PCA is used on GMV data for dimensionality reduction.
#
# Key differences between tasks:
#  1) Math uses residuals of "numop_std" and "mathprob_std"
#     (after regressing out SES).
#  2) Reading uses residuals of "wordread_std" and "readcomp_std"
#     (after regressing out SES).
#
# The script loops over both tasks.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear environment
rm(list = ls())

# Load custom utility functions
source("scripts/utility/myRFunc.R")

# ------------------------------------------------------------------
# 2. File Paths and PCA Settings
# ------------------------------------------------------------------

# Input data files
bn_atlas_file <- "data/atlas/bn_atlas.xlsx"
beh_file      <- "data/subjectlist/subjectlist_cmi_n760.csv"
gmv_file      <- "data/gmv/gmv_cmi_n760.csv"

# PCA configuration
apply_pca <- 1   # 1 = apply PCA to brain data
prefix    <- "CCA_PCA"
condition <- "roi_gmv_brainnetome"

# ------------------------------------------------------------------
# 3. Define Analyses (Reading and Math)
# ------------------------------------------------------------------

analyses <- list(
  math = list(
    output_path = "results/cca/cmi/wholebrain_cca_cmi_math_ctrSES",
    postfix     = "mathstd_ageinmodel",
    cols        = c("age", "numop_std", "mathprob_std", "Barratt_Total")
  ),
  reading = list(
    output_path = "results/cca/cmi/wholebrain_cca_cmi_reading_ctrSES",
    postfix     = "readstd_ageinmodel",
    cols        = c("age", "wordread_std", "readcomp_std", "Barratt_Total")
  )
)

# ------------------------------------------------------------------
# 4. Data Loading and Preprocessing
# ------------------------------------------------------------------

# Load data
data <- load_data(bn_atlas_file, beh_file, gmv_file)

# Verify participant IDs match
sum(data$beh$oak_id == data$pids)

# --------------------
# Outlier Detection
# --------------------

beh_df <- as.data.frame(data$beh)
ses <- beh_df$Barratt_Total  # SES variable

# Identify behavioral outliers (age, math, reading)
col_of_interest <- c("age", "mathprob_std", "numop_std", "readcomp_std", "wordread_std")
beh_df <- beh_df[, which(colnames(beh_df) %in% col_of_interest)]
beh_outliers <- lapply(beh_df, identify_outliers)
beh_exo_idx  <- sort(unique(unlist(beh_outliers)))

# Identify brain outliers
brain_outliers <- lapply(as.data.frame(data$brain), identify_outliers)
brain_exo_idx  <- sort(unique(unlist(brain_outliers)))

# Identify NA values in SES
na_idx <- which(is.na(ses))

# Combine all exclusions
exo_idx <- unique(c(beh_exo_idx, brain_exo_idx, na_idx))

# Remove outliers and NAs
if (length(exo_idx) > 0) {
  data$beh   <- data$beh[-exo_idx, ]
  data$brain <- data$brain[-exo_idx, ]
  data$pids  <- data$pids[-exo_idx]
  data$numsub <- data$numsub - length(exo_idx)
}

# --------------------
# PCA on Brain Data
# --------------------
if (apply_pca == 1) {
  my.pca <- mypca(data$brain, 0.9)  # Keep 90% variance
}
Xorig <- data$brain

# ------------------------------------------------------------------
# 5. Loop Over Analyses (Reading and Math)
# ------------------------------------------------------------------

for (task in names(analyses)) {
  cat("\n----------------------------\n")
  cat("Running CCA for:", task, " (SES-Controlled)\n")
  cat("----------------------------\n")
  
  # Extract parameters
  output_path <- analyses[[task]]$output_path
  postfix     <- analyses[[task]]$postfix
  cols        <- analyses[[task]]$cols
  
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare brain matrix (PCA if needed)
  if (apply_pca == 1) {
    X <- my.pca$x[, 1:my.pca$numPC]
  } else {
    X <- data$brain
  }
  
  # Extract behavioral data
  YY <- data$beh[, which(colnames(data$beh) %in% cols)]
  
  # ----------------------------------------
  # Regress out SES (Barratt_Total)
  # ----------------------------------------
  if (task == "math") {
    numop_resid    <- resid(lm(numop_std ~ Barratt_Total, data = YY))
    mathprob_resid <- resid(lm(mathprob_std ~ Barratt_Total, data = YY))
    Y <- data.frame(
      age = YY$age,
      numop_std_resid = numop_resid,
      mathprob_std_resid = mathprob_resid
    )
  } else if (task == "reading") {
    wordread_resid <- resid(lm(wordread_std ~ Barratt_Total, data = YY))
    readcomp_resid <- resid(lm(readcomp_std ~ Barratt_Total, data = YY))
    Y <- data.frame(
      age = YY$age,
      wordread_std_resid = wordread_resid,
      readcomp_std_resid = readcomp_resid
    )
  }
  
  # -----------------------
  # Run CCA with 5000 permutations
  # -----------------------
  res <- mycca_noscale(X, Y, 5000)
  
  # -----------------------
  # Print Key CCA Results
  # -----------------------
  cat("\n### Multivariate Test (Pillai's Trace) ###\n")
  cat("Pillai's trace:", res$Pillai, "\n")
  cat("Pillai's p-value:", res$pillai.p, "\n")
  cat("Pillai's p-value (permutation):", res$pillai.p.perm, "\n")
  
  cat("\n### Canonical Correlations ###\n")
  cat("Canonical correlations:", paste(res$cancor, collapse = ", "), "\n")
  
  cat("\n### Significance of Canonical Correlations ###\n")
  cat("Permutation p-values:", paste(res$cancor.p.perm, collapse = ", "), "\n")
  cat("Adjusted permutation p-values:", paste(res$cancor.p.perm.adjust, collapse = ", "), "\n")
  
  cat("\n### Wilks' Lambda Dimension Test ###\n")
  p.asym(res$cancor, nrow(X), ncol(X), ncol(Y), tstat = "Wilks")
  
  # -----------------------
  # Save Results
  # -----------------------
  if (apply_pca == 1) {
    fname <- paste(output_path, '/', prefix, '_', condition, '_', postfix,
                   '_perm5000_pcaNoscale_ccaNoScale.RData', sep = "")
    save(data, Xorig, X, Y, res, my.pca, file = fname)
  } else {
    fname <- paste(output_path, '/', prefix, '_', condition, '_', postfix,
                   '_perm5000_ccaNoScale.RData', sep = "")
    save(data, Xorig, X, Y, res, file = fname)
  }
  
  # -----------------------
  # Save Weights and Plots
  # -----------------------
  write_weights(res, Xorig, Y, apply_pca, output_path, prefix, condition, postfix)
  myplot_YCy(res, output_path, prefix, condition, postfix)
  myplot_YcoefStd(res, output_path, prefix, condition, postfix)
  myplot_cca_scatter(res, output_path, prefix, condition, postfix)
  myplot_XcoefStd(res, output_path, prefix, condition, postfix)
  myplot_XCx(res, output_path, prefix, condition, postfix)
}

cat("\nAll analyses (Reading and Math) with SES control completed!\n")
