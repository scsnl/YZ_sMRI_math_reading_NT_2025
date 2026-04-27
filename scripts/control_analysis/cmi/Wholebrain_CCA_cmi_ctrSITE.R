############################################################
# Canonical Correlation Analysis (CCA) - CMI Cohort
# Author: Yuan Zhang
# Date: 2026-03-24
#
# Description:
# This script runs CCA analyses for both math and reading
# tasks on the CMI cohort, controlling for site on the
# brain side only. Specifically, SITE is regressed out
# from ROI GMV data before CCA. PCA is used on GMV data
# for dimensionality reduction.
#
# Key differences between tasks:
#  1) Math uses age, "numop_std", and "mathprob_std"
#  2) Reading uses age, "wordread_std", and "readcomp_std"
#
# The script loops over both tasks.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory (change to your local project folder)
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear environment
rm(list = ls())

# Load custom utility functions: load_data, identify_outliers, mypca, mycca_noscale, etc.
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

# Parameters for each task
analyses <- list(
  math = list(
    output_path = "results/cca/cmi/wholebrain_cca_cmi_math_ctrSite",
    postfix     = "mathstd_ageinmodel",
    cols        = c("age", "numop_std", "mathprob_std")
  ),
  reading = list(
    output_path = "results/cca/cmi/wholebrain_cca_cmi_reading_ctrSite",
    postfix     = "readstd_ageinmodel",
    cols        = c("age", "wordread_std", "readcomp_std")
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

# Identify outliers in brain GMV data
brain_outliers <- lapply(as.data.frame(data$brain), identify_outliers)
brain_exo_idx  <- sort(unique(unlist(brain_outliers)))

# Identify outliers in behavioral data
beh_df <- as.data.frame(data$beh)
col_of_interest <- c("age", "mathprob_std", "numop_std", "readcomp_std", "wordread_std")
beh_df <- beh_df[, which(colnames(beh_df) %in% col_of_interest)]
beh_outliers <- lapply(beh_df, identify_outliers)
beh_exo_idx  <- sort(unique(unlist(beh_outliers)))

# Combine outliers and remove them
exo_idx <- unique(c(beh_exo_idx, brain_exo_idx))
data$beh   <- data$beh[-exo_idx, ]
data$brain <- data$brain[-exo_idx, ]
data$pids  <- data$pids[-exo_idx]
data$numsub <- data$numsub - length(exo_idx)

# --------------------
# Regress out SITE from brain data only
# --------------------
# SITE should be a categorical variable in data$beh
data$beh$SITE <- as.factor(data$beh$SITE)

brain_site_resid <- sapply(1:ncol(data$brain), function(i) {
  resid(lm(data$brain[, i] ~ data$beh$SITE))
})

brain_site_resid <- as.matrix(brain_site_resid)
colnames(brain_site_resid) <- colnames(data$brain)
rownames(brain_site_resid) <- rownames(data$brain)

# Keep original brain data for saving / plotting weights
Xorig <- data$brain

# --------------------
# PCA on Site-adjusted Brain Data
# --------------------
if (apply_pca == 1) {
  my.pca <- mypca(brain_site_resid, 0.9)  # Keep 90% variance
}

# ------------------------------------------------------------------
# 5. Loop Over Analyses (Reading and Math)
# ------------------------------------------------------------------

for (task in names(analyses)) {
  cat("\n----------------------------\n")
  cat("Running CCA for:", task, "\n")
  cat("----------------------------\n")
  
  # Extract parameters
  output_path <- analyses[[task]]$output_path
  postfix     <- analyses[[task]]$postfix
  cols        <- analyses[[task]]$cols
  
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare data for CCA
  if (apply_pca == 1) {
    X <- my.pca$x[, 1:my.pca$numPC]
  } else {
    X <- brain_site_resid
  }
  
  Y <- data$beh[, which(colnames(data$beh) %in% cols)]
  
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
    save(data, Xorig, X, Y, res, my.pca, brain_site_resid, file = fname)
  } else {
    fname <- paste(output_path, '/', prefix, '_', condition, '_', postfix,
                   '_perm5000_ccaNoScale.RData', sep = "")
    save(data, Xorig, X, Y, res, brain_site_resid, file = fname)
  }
  
  # -----------------------
  # Save Weights and Plots
  # -----------------------
  # Note: write_weights uses Xorig, so weight maps remain in original ROI space
  write_weights(res, Xorig, Y, apply_pca, output_path, prefix, condition, postfix)
  myplot_YCy(res, output_path, prefix, condition, postfix)
  myplot_YcoefStd(res, output_path, prefix, condition, postfix)
  myplot_cca_scatter(res, output_path, prefix, condition, postfix)
  myplot_XcoefStd(res, output_path, prefix, condition, postfix)
  myplot_XCx(res, output_path, prefix, condition, postfix)
}

cat("\nAll analyses (Reading and Math) with site control completed!\n")