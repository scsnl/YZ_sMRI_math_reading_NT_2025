############################################################
# Compare GMV weight maps across models within each cohort
#
# For each cohort, compare the Mode 2 GMV weight maps from:
#   1) math-alone model # mode2 = math mode
#   2) reading-alone model # mode2 = reading mode
#   3) math+reading combined model # mode2 = shared mode
#
# GMV weights are computed as:
#   cor(Xorig, -1 * res$Cx[, 2])
#
# Author: Yuan Zhang
# Date: 2026-04-03
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(Hmisc)
library(corrplot)

# ------------------------------------------------------------------
# 2. Helper function: load one model and compute Mode 2 GMV weights
# ------------------------------------------------------------------
get_mode2_gmv_weights <- function(rdata_file, mode_idx = 2, flip_sign = TRUE) {
  load(rdata_file)
  
  # GMV weight map from original ROI-space brain data
  w <- as.vector(cor(Xorig, res$Cx[, mode_idx]))
  
  if (flip_sign) {
    w <- -1 * w
  }
  
  return(w)
}

# ------------------------------------------------------------------
# 3. Helper function: compare maps within one cohort
# ------------------------------------------------------------------
compare_within_cohort <- function(math_file, read_file, combined_file,
                                  cohort_name, out_dir,
                                  mode_math = 2, mode_read = 2, mode_combined = 2,
                                  flip_sign = TRUE) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load weight maps
  w_math <- get_mode2_gmv_weights(math_file, mode_idx = mode_math, flip_sign = flip_sign)
  w_read <- get_mode2_gmv_weights(read_file, mode_idx = mode_read, flip_sign = flip_sign)
  w_comb <- get_mode2_gmv_weights(combined_file, mode_idx = mode_combined, flip_sign = flip_sign)
  
  # Put into one dataframe
  df <- data.frame(
    Math_Alone = w_math,
    Reading_Alone = w_read,
    Combined_Mode2 = w_comb
  )
  
  # Save raw weights
  write.csv(
    df,
    file.path(out_dir, paste0(cohort_name, "_mode2_gmv_weights.csv")),
    row.names = FALSE
  )
  
  # Correlation matrix and p-value matrix
  rc <- rcorr(as.matrix(df), type = "pearson")
  cor_mat <- rc$r
  p_mat <- rc$P
  
  write.csv(
    cor_mat,
    file.path(out_dir, paste0(cohort_name, "_mode2_gmv_weight_correlations.csv"))
  )
  
  write.csv(
    p_mat,
    file.path(out_dir, paste0(cohort_name, "_mode2_gmv_weight_pvalues.csv"))
  )
  
  # Pairwise tests
  pairs <- combn(colnames(df), 2, simplify = FALSE)
  pairwise_results <- do.call(rbind, lapply(pairs, function(pair) {
    test <- cor.test(df[[pair[1]]], df[[pair[2]]], method = "pearson")
    data.frame(
      cohort = cohort_name,
      var1 = pair[1],
      var2 = pair[2],
      r = unname(test$estimate),
      p = test$p.value,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(
    pairwise_results,
    file.path(out_dir, paste0(cohort_name, "_mode2_gmv_weight_pairwise_tests.csv")),
    row.names = FALSE
  )
  
  # Heatmap
  eps_file <- file.path(out_dir, paste0(cohort_name, "_mode2_gmv_weight_heatmap.eps"))
  postscript(eps_file, width = 5.2, height = 5.2, horizontal = FALSE,
             onefile = FALSE, paper = "special")
  corrplot(
    cor_mat,
    method = "color",
    col = colorRampPalette(c("orange", "tomato"))(200),
    type = "upper",
    diag = FALSE,
    addCoef.col = "white",
    number.cex =  0.85,
    tl.col = "black",
    tl.srt = 45,
    cl.pos = "n",
    is.corr = FALSE,
    cl.lim = c(0, 1),
    mar = c(0, 0, 1, 0)
  )
  title(main = paste0(cohort_name, ": Correlation of Mode 2 GMV Weight Maps"))
  dev.off()
  
  cat("\n====================================\n")
  cat("Cohort:", cohort_name, "\n")
  cat("====================================\n")
  print(round(cor_mat, 4))
  cat("\nPairwise tests:\n")
  print(pairwise_results)
  
  invisible(list(
    weights = df,
    cor_mat = cor_mat,
    p_mat = p_mat,
    pairwise = pairwise_results
  ))
}

# ------------------------------------------------------------------
# 4. File paths
# ------------------------------------------------------------------

# ---- CMI ----
cmi_math_file <- paste0(
  "results/cca/cmi/wholebrain_cca_cmi_math/",
  "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)

cmi_read_file <- paste0(
  "results/cca/cmi/wholebrain_cca_cmi_reading/",
  "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)

cmi_combined_file <- paste0(
  "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/",
  "CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)

# ---- Stanford ----
stan_math_file <- paste0(
  "results/cca/stanford/wholebrain_cca_stanford_math/",
  "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)

stan_read_file <- paste0(
  "results/cca/stanford/wholebrain_cca_stanford_reading/",
  "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)

stan_combined_file <- paste0(
  "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined/",
  "CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)

# ------------------------------------------------------------------
# 5. Run comparisons
# ------------------------------------------------------------------

# Set flip_sign = TRUE if you want all maps sign-flipped before comparing
res_cmi <- compare_within_cohort(
  math_file = cmi_math_file,
  read_file = cmi_read_file,
  combined_file = cmi_combined_file,
  cohort_name = "CMI",
  out_dir = "results/cca/cmi/model_weight_map_comparison/",
  mode_math = 2,
  mode_read = 2,
  mode_combined = 2,
  flip_sign = TRUE
)

res_stan <- compare_within_cohort(
  math_file = stan_math_file,
  read_file = stan_read_file,
  combined_file = stan_combined_file,
  cohort_name = "Stanford",
  out_dir = "results/cca/stanford/model_weight_map_comparison/",
  mode_math = 2,
  mode_read = 2,
  mode_combined = 2,
  flip_sign = TRUE
)