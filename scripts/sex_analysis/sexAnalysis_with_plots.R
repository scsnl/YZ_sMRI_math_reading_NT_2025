##########################################################################
# Sex Analysis of CCA GMV Weights for CMI-HBN and Stanford Cohorts
# Author: Yuan Zhang
# Date: 2026-03-24
#
# Description:
#   1. Loads math and reading CCA results for both CMI-HBN and Stanford.
#   2. Extracts the brain-side CCA scores (U) and behavior-side CCA scores (V)
#      for math and reading.
#   3. Computes GMV weight maps (full sample and both sex).
#   4. Saves GMV weight maps to CSV files for both math and reading.
#   5. Creates correlation plots comparing GMV weight maps across sex.
#   6. Computes canonical correlation for each sex.
#   7. Examines whether there are sex differences in canonical correlations.
#   8. Creates scatter plots of U-V association by sex.
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
    output_dir = "results/sex_analysis/cmi/"
  ),
  Stanford = list(
    math_path = "results/cca/stanford/wholebrain_cca_stanford_math/",
    math_file = "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    read_path = "results/cca/stanford/wholebrain_cca_stanford_reading/",
    read_file = "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_dir = "results/sex_analysis/stanford/"
  )
)

# Ensure output directories exist
for (ds in names(datasets)) {
  dir.create(datasets[[ds]]$output_dir, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------------
# 3. Helper Functions
# ------------------------------------------------------------------

# Fisher z test for difference between two independent correlations
compare_correlations_fisher <- function(r1, n1, r2, n2) {
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z  <- (z1 - z2) / se
  p  <- 2 * (1 - pnorm(abs(z)))
  data.frame(
    r1 = r1, n1 = n1,
    r2 = r2, n2 = n2,
    z = z, p = p
  )
}

# Extract sex labels consistently across datasets
get_sex_vector <- function(beh, dataset_name) {
  if (dataset_name == "CMI") {
    sex <- ifelse(beh$Sex == "M", "Male", "Female")
  } else {
    sex <- ifelse(beh$genderF == 0, "Male", "Female")
  }
  return(sex)
}

# Pairwise cor.test results for columns in a dataframe
get_pairwise_cor_tests <- function(df, analysis_name) {
  pairs <- combn(colnames(df), 2, simplify = FALSE)
  out <- lapply(pairs, function(pair) {
    test <- cor.test(df[[pair[1]]], df[[pair[2]]], method = "pearson")
    data.frame(
      analysis = analysis_name,
      var1 = pair[1],
      var2 = pair[2],
      r = unname(test$estimate),
      p = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# Corrplot and save
plot_corr <- function(df, output_file) {
  stats <- rcorr(as.matrix(df))  # Pearson correlation matrix
  cor_matrix <- stats$r
  colnames(cor_matrix) <- rownames(cor_matrix) <- colnames(df)
  
  my_col <- colorRampPalette(c("orange", "tomato"))(200)
  
  corrplot(
    cor_matrix, method = "color", col = my_col, type = "upper", diag = FALSE,
    addCoef.col = "white", number.cex = 0.85,
    tl.col = "black", tl.srt = 45, cl.pos = "n", # col.lim = c(0.75, 1),
    is.corr = FALSE, mar = c(0, 0, 1, 0)
  )
  dev.copy2eps(file = output_file, width = 5, height = 5)
  dev.off()
}

# Scatter plot for U-V association by sex
plot_uv_by_sex <- function(U, V, sex, title_text, output_file) {
  df <- data.frame(U = U, V = V, Sex = sex)
  
  # Overall and group-specific statistics
  cor_m <- cor.test(df$U[df$Sex == "Male"], df$V[df$Sex == "Male"])
  cor_f <- cor.test(df$U[df$Sex == "Female"], df$V[df$Sex == "Female"])
  
  label_text <- paste0(
    "Male: r = ", sprintf("%.3f", unname(cor_m$estimate)),
    ", p = ", formatC(cor_m$p.value, format = "e", digits = 2), "\n",
    "Female: r = ", sprintf("%.3f", unname(cor_f$estimate)),
    ", p = ", formatC(cor_f$p.value, format = "e", digits = 2)
  )
  
  p <- ggplot(df, aes(x = U, y = V, color = Sex)) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    annotate(
      "text",
      x = min(df$U, na.rm = TRUE),
      y = max(df$V, na.rm = TRUE),
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = title_text,
      x = "Brain canonical variate (U)",
      y = "Behavior canonical variate (V)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  ggsave(output_file, p, width = 6, height = 5)
}

# ------------------------------------------------------------------
# 4. Function to Run Sex Analysis
# ------------------------------------------------------------------
analyze_sex <- function(math_path, math_file, read_path, read_file, output_dir, dataset_name) {
  
  cat("\n===================================\n")
  cat("Processing Dataset:", dataset_name, "\n")
  cat("===================================\n")
  
  # --- Load Math CCA Results ---
  load(file.path(math_path, math_file))
  Umath <- res$Cx[, 2]   # Brain-side CCA mode 2 scores
  Vmath <- res$Cy[, 2]   # Behavior-side CCA mode 2 scores
  gmv   <- data$brain    # GMV data (nsub x nroi)
  beh   <- data$beh      # Behavioral data
  
  # Keep math data references before loading reading
  gmv_math <- gmv
  beh_math <- beh
  
  # --- Load Reading CCA Results ---
  load(file.path(read_path, read_file))
  Uread <- res$Cx[, 2]   # Brain-side CCA mode 2 scores
  Vread <- res$Cy[, 2]   # Behavior-side CCA mode 2 scores
  gmv_read <- data$brain
  beh_read <- data$beh
  
  # ------------------------------------------------------------------
  # 5. Compute GMV Weight Maps (Full Sample & Each Sex)
  # ------------------------------------------------------------------
  xloading_math_full <- as.vector(cor(gmv_math, -1 * Umath))
  xloading_read_full <- as.vector(cor(gmv_read, -1 * Uread))
  
  sex_math <- get_sex_vector(beh_math, dataset_name)
  sex_read <- get_sex_vector(beh_read, dataset_name)
  
  idx_m_math <- which(sex_math == "Male")
  idx_f_math <- which(sex_math == "Female")
  idx_m_read <- which(sex_read == "Male")
  idx_f_read <- which(sex_read == "Female")
  
  xloading_math_m <- as.vector(cor(gmv_math[idx_m_math, ], -1 * Umath[idx_m_math]))
  xloading_math_f <- as.vector(cor(gmv_math[idx_f_math, ], -1 * Umath[idx_f_math]))
  
  xloading_read_m <- as.vector(cor(gmv_read[idx_m_read, ], -1 * Uread[idx_m_read]))
  xloading_read_f <- as.vector(cor(gmv_read[idx_f_read, ], -1 * Uread[idx_f_read]))
  
  df_math <- data.frame(
    Full = xloading_math_full,
    Male = xloading_math_m,
    Female = xloading_math_f
  )
  
  df_read <- data.frame(
    Full = xloading_read_full,
    Male = xloading_read_m,
    Female = xloading_read_f
  )
  
  write.csv(df_math, file.path(output_dir, "gmv_weights_sex_math.csv"), row.names = FALSE)
  write.csv(df_read, file.path(output_dir, "gmv_weights_sex_read.csv"), row.names = FALSE)
  
  # ------------------------------------------------------------------
  # 6. Correlations between GMV maps + p-values
  # ------------------------------------------------------------------
  cor_tests_math <- get_pairwise_cor_tests(df_math, analysis_name = "Math_GMV_Weights")
  cor_tests_read <- get_pairwise_cor_tests(df_read, analysis_name = "Reading_GMV_Weights")
  cor_tests_all  <- rbind(cor_tests_math, cor_tests_read)
  
  write.csv(cor_tests_all, file.path(output_dir, "gmv_weights_map_correlations_with_p.csv"), row.names = FALSE)
  
  # ------------------------------------------------------------------
  # 7. Correlation Plots for Math & Reading GMV Weights
  # ------------------------------------------------------------------
  plot_corr(df_math, file.path(output_dir, "gmv_weights_correlation_math.eps"))
  plot_corr(df_read, file.path(output_dir, "gmv_weights_correlation_read.eps"))
  
  # ------------------------------------------------------------------
  # 8. Compute canonical correlation for each sex
  # ------------------------------------------------------------------
  cor_math_m <- cor.test(Umath[idx_m_math], Vmath[idx_m_math], method = "pearson")
  cor_math_f <- cor.test(Umath[idx_f_math], Vmath[idx_f_math], method = "pearson")
  cor_read_m <- cor.test(Uread[idx_m_read], Vread[idx_m_read], method = "pearson")
  cor_read_f <- cor.test(Uread[idx_f_read], Vread[idx_f_read], method = "pearson")
  
  canon_corr_df <- data.frame(
    analysis = c("Math", "Math", "Reading", "Reading"),
    sex = c("Male", "Female", "Male", "Female"),
    n = c(length(idx_m_math), length(idx_f_math), length(idx_m_read), length(idx_f_read)),
    r = c(unname(cor_math_m$estimate), unname(cor_math_f$estimate),
          unname(cor_read_m$estimate), unname(cor_read_f$estimate)),
    p = c(cor_math_m$p.value, cor_math_f$p.value,
          cor_read_m$p.value, cor_read_f$p.value)
  )
  
  write.csv(canon_corr_df, file.path(output_dir, "canonical_correlations_by_sex.csv"), row.names = FALSE)
  
  # ------------------------------------------------------------------
  # 9. Examine sex differences in canonical correlations
  # ------------------------------------------------------------------
  math_sex_diff <- compare_correlations_fisher(
    r1 = unname(cor_math_m$estimate), n1 = length(idx_m_math),
    r2 = unname(cor_math_f$estimate), n2 = length(idx_f_math)
  )
  math_sex_diff$analysis <- "Math"
  
  read_sex_diff <- compare_correlations_fisher(
    r1 = unname(cor_read_m$estimate), n1 = length(idx_m_read),
    r2 = unname(cor_read_f$estimate), n2 = length(idx_f_read)
  )
  read_sex_diff$analysis <- "Reading"
  
  sex_diff_df <- rbind(math_sex_diff, read_sex_diff)
  sex_diff_df <- sex_diff_df[, c("analysis", "r1", "n1", "r2", "n2", "z", "p")]
  colnames(sex_diff_df) <- c("analysis", "r_male", "n_male", "r_female", "n_female", "z", "p")
  
  write.csv(sex_diff_df, file.path(output_dir, "sex_difference_in_canonical_correlations.csv"), row.names = FALSE)
  
  # ------------------------------------------------------------------
  # 10. Scatter plots by sex
  # ------------------------------------------------------------------
  plot_uv_by_sex(
    U = Umath, V = Vmath, sex = sex_math,
    title_text = paste0(dataset_name, ": Math CCA U-V Association by Sex"),
    output_file = file.path(output_dir, "scatter_math_uv_by_sex.png")
  )
  
  plot_uv_by_sex(
    U = Uread, V = Vread, sex = sex_read,
    title_text = paste0(dataset_name, ": Reading CCA U-V Association by Sex"),
    output_file = file.path(output_dir, "scatter_read_uv_by_sex.png")
  )
  
  # ------------------------------------------------------------------
  # 11. Save a simple text summary
  # ------------------------------------------------------------------
  sink(file.path(output_dir, "sex_analysis_summary.txt"))
  
  cat("Dataset:", dataset_name, "\n\n")
  
  cat("=== GMV weight map correlations ===\n")
  print(cor_tests_all)
  cat("\n")
  
  cat("=== Canonical correlations by sex ===\n")
  print(canon_corr_df)
  cat("\n")
  
  cat("=== Sex differences in canonical correlations (Fisher z test) ===\n")
  print(sex_diff_df)
  cat("\n")
  
  sink()
  
  cat("Finished:", dataset_name, "\n")
}

# ------------------------------------------------------------------
# 12. Run Analysis for Both Datasets
# ------------------------------------------------------------------
for (ds in names(datasets)) {
  analyze_sex(
    math_path = datasets[[ds]]$math_path,
    math_file = datasets[[ds]]$math_file,
    read_path = datasets[[ds]]$read_path,
    read_file = datasets[[ds]]$read_file,
    output_dir = datasets[[ds]]$output_dir,
    dataset_name = ds
  )
}