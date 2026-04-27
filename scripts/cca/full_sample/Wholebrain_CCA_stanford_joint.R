############################################################
# Canonical Correlation Analysis (CCA) - Stanford Cohort
# Author: Yuan Zhang
# Date: 2026-04-02
#
# Description:
# This script runs a combined CCA analysis for the Stanford
# cohort, including both math- and reading-related behavioral
# measures, while controlling for age. PCA is used on GMV
# data for dimensionality reduction.
#
# Behavioral variables included:
#   1) age
#   2) numop_std
#   3) mathrea_std
#   4) wordread_std
#   5) readcomp_std
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
beh_file      <- "data/subjectlist/subjectlist_stanford_n231.csv"
gmv_file      <- "data/gmv/gmv_stanford_n231.csv"

# PCA configuration
apply_pca <- 1   # 1 = apply PCA to brain data
prefix    <- "CCA_PCA"
condition <- "roi_gmv_brainnetome"

# Output settings
output_path <- "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined"
postfix     <- "mathreadingstd_ageinmodel"
cols        <- c("age", "numop_std", "mathrea_std", "wordread_std", "readcomp_std")

# ------------------------------------------------------------------
# 3. Data Loading and Preprocessing
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
col_of_interest <- c("age", "mathrea_std", "numop_std", "readcomp_std", "wordread_std")
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
# PCA on Brain Data
# --------------------
if (apply_pca == 1) {
  my.pca <- mypca(data$brain, 0.9)  # Keep 90% variance
}
Xorig <- data$brain

# ------------------------------------------------------------------
# 4. Prepare Data for CCA
# ------------------------------------------------------------------

# Create output directory
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# Prepare brain matrix
if (apply_pca == 1) {
  X <- my.pca$x[, 1:my.pca$numPC]
} else {
  X <- data$brain
}

# Prepare behavior matrix
# Y <- data$beh[, which(colnames(data$beh) %in% cols)]
Y <- data$beh[, cols, drop = FALSE]

# ------------------------------------------------------------------
# 5. Run CCA with 5000 permutations
# ------------------------------------------------------------------
cat("\n----------------------------\n")
cat("Running combined CCA for math + reading\n")
cat("----------------------------\n")

res <- mycca_noscale(X, Y, 5000)

# ------------------------------------------------------------------
# 6. Print Key CCA Results
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# 7. Save Results
# ------------------------------------------------------------------
if (apply_pca == 1) {
  fname <- paste(output_path, '/', prefix, '_', condition, '_', postfix,
                 '_perm5000_pcaNoscale_ccaNoScale.RData', sep = "")
  save(data, Xorig, X, Y, res, my.pca, file = fname)
} else {
  fname <- paste(output_path, '/', prefix, '_', condition, '_', postfix,
                 '_perm5000_ccaNoScale.RData', sep = "")
  save(data, Xorig, X, Y, res, file = fname)
}

# ------------------------------------------------------------------
# 8. Save Weights and Plots
# ------------------------------------------------------------------
write_weights(res, Xorig, Y, apply_pca, output_path, prefix, condition, postfix)
myplot_YCy(res, output_path, prefix, condition, postfix)
myplot_YcoefStd(res, output_path, prefix, condition, postfix)
myplot_cca_scatter(res, output_path, prefix, condition, postfix)
myplot_XcoefStd(res, output_path, prefix, condition, postfix)
myplot_XCx(res, output_path, prefix, condition, postfix)

cat("\nCombined CCA analysis (Math + Reading) completed!\n")