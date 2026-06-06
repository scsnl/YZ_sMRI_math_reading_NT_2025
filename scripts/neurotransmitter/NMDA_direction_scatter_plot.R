############################################################
# NMDA scatter plots for CCA-derived GMV maps
# Author: Yuan Zhang
# Date: 2026-05-29
#
# Description:
# This script plots associations between CCA-derived GMV
# structural phenotypes and the NMDA receptor map for:
#   1. CMI math
#   2. CMI reading
#   3. Stanford math
#   4. Stanford reading
#
# Saves:
#   - scatter plots: EPS / PDF / PNG
#   - source data for scatter points
#   - source data for fitted regression line and 95% CI band
#   - Pearson r summary table
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ggplot2)

out_dir <- "results/neurotransmitter/NMDA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------

format_r <- function(r) {
  sprintf("%.2f", r)
}

# Extract CCA GMV map from original ROI-space GMV data.
# This is the correlation between each ROI's original GMV and the brain canonical variate.
get_mode_gmv_weights <- function(rdata_file, mode_idx = 2, flip_sign = TRUE) {
  load(rdata_file)
  
  w <- as.vector(cor(Xorig, res$Cx[, mode_idx]))
  
  if (flip_sign) {
    w <- -1 * w
  }
  
  return(w)
}

plot_nmda_scatter <- function(
    rdata_file,
    output_prefix,
    plot_title,
    x_label,
    mode_idx = 2,
    flip_sign = TRUE
) {
  
  # -----------------------------
  # Load GMV structural phenotype
  # -----------------------------
  gmv_weight <- get_mode_gmv_weights(
    rdata_file = rdata_file,
    mode_idx = mode_idx,
    flip_sign = flip_sign
  )
  
  # -----------------------------
  # Load NMDA receptor map
  # -----------------------------
  fname_receptors <- "data/neurotransmitter/receptor_data_bn246.csv"
  receptors <- read.csv(fname_receptors, header = FALSE)
  
  # Receptor order:
  # c('5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT',
  #   'A4B2', 'CB1', 'D1', 'D2', 'DAT', 'GABAa', 'H3',
  #   'M1', 'mGluR5', 'MOR', 'NET', 'NMDA', 'VAChT')
  nmda_map <- receptors[1:218, 18]
  
  df_scatter <- data.frame(
    ROI = seq_len(218),
    GMV_structural_phenotype = as.numeric(gmv_weight),
    NMDA = as.numeric(nmda_map)
  )
  
  df_scatter <- df_scatter[complete.cases(df_scatter), ]
  
  # -----------------------------
  # Pearson r for visual annotation
  # -----------------------------
  r_value <- cor(
    df_scatter$GMV_structural_phenotype,
    df_scatter$NMDA,
    method = "pearson"
  )
  
  annotation_text <- paste0("r = ", format_r(r_value))
  
  # -----------------------------
  # Save scatter source data
  # -----------------------------
  write.csv(
    df_scatter,
    file.path(out_dir, paste0(output_prefix, "_source_data_scatter_points.csv")),
    row.names = FALSE
  )
  
  # -----------------------------
  # Fitted regression line and 95% CI band
  # -----------------------------
  lm_fit <- lm(NMDA ~ GMV_structural_phenotype, data = df_scatter)
  
  x_grid <- data.frame(
    GMV_structural_phenotype = seq(
      min(df_scatter$GMV_structural_phenotype, na.rm = TRUE),
      max(df_scatter$GMV_structural_phenotype, na.rm = TRUE),
      length.out = 200
    )
  )
  
  pred_ci <- predict(
    lm_fit,
    newdata = x_grid,
    interval = "confidence",
    level = 0.95
  )
  
  df_fit_ci <- data.frame(
    GMV_structural_phenotype = x_grid$GMV_structural_phenotype,
    fitted_NMDA = pred_ci[, "fit"],
    ci_lower = pred_ci[, "lwr"],
    ci_upper = pred_ci[, "upr"]
  )
  
  write.csv(
    df_fit_ci,
    file.path(out_dir, paste0(output_prefix, "_source_data_fit_95CI.csv")),
    row.names = FALSE
  )
  
  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(
    df_scatter,
    aes(x = GMV_structural_phenotype, y = NMDA)
  ) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = lm, se = TRUE) +
    annotate(
      "text",
      x = min(df_scatter$GMV_structural_phenotype, na.rm = TRUE),
      y = max(df_scatter$NMDA, na.rm = TRUE),
      label = annotation_text,
      hjust = 0,
      vjust = 1,
      size = 5
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = plot_title,
      x = x_label,
      y = "NMDA receptor map"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save plots
  ggsave(
    filename = file.path(out_dir, paste0(output_prefix, ".eps")),
    plot = p,
    units = "in",
    width = 5,
    height = 5,
    device = cairo_ps
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(output_prefix, ".pdf")),
    plot = p,
    units = "in",
    width = 5,
    height = 5
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(output_prefix, ".png")),
    plot = p,
    units = "in",
    width = 5,
    height = 5,
    dpi = 300
  )
  
  # Return summary
  summary_df <- data.frame(
    analysis = output_prefix,
    mode = mode_idx,
    n_ROI = nrow(df_scatter),
    pearson_r = r_value,
    source_data_scatter_points = file.path(
      out_dir,
      paste0(output_prefix, "_source_data_scatter_points.csv")
    ),
    source_data_fit_95CI = file.path(
      out_dir,
      paste0(output_prefix, "_source_data_fit_95CI.csv")
    ),
    stringsAsFactors = FALSE
  )
  
  return(list(plot = p, scatter_data = df_scatter, fit_data = df_fit_ci, summary = summary_df))
}

# ------------------------------------------------------------------
# 3. Define analyses
# ------------------------------------------------------------------

analyses <- list(
  cmi_math = list(
    rdata_file = "results/cca/cmi/wholebrain_cca_cmi_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_prefix = "scatter_cmi_math_NMDA",
    plot_title = "CMI-HBN math",
    x_label = "Math-related GMV structural phenotype"
  ),
  
  cmi_reading = list(
    rdata_file = "results/cca/cmi/wholebrain_cca_cmi_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_prefix = "scatter_cmi_reading_NMDA",
    plot_title = "CMI-HBN reading",
    x_label = "Reading-related GMV structural phenotype"
  ),
  
  stanford_math = list(
    rdata_file = "results/cca/stanford/wholebrain_cca_stanford_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_prefix = "scatter_stanford_math_NMDA",
    plot_title = "Stanford math",
    x_label = "Math-related GMV structural phenotype"
  ),
  
  stanford_reading = list(
    rdata_file = "results/cca/stanford/wholebrain_cca_stanford_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    output_prefix = "scatter_stanford_reading_NMDA",
    plot_title = "Stanford reading",
    x_label = "Reading-related GMV structural phenotype"
  )
)

# ------------------------------------------------------------------
# 4. Run analyses
# ------------------------------------------------------------------

all_results <- list()

for (name in names(analyses)) {
  cat("\n----------------------------------------\n")
  cat("Running:", name, "\n")
  cat("----------------------------------------\n")
  
  a <- analyses[[name]]
  
  if (!file.exists(a$rdata_file)) {
    warning(paste("File not found, skipping:", a$rdata_file))
    next
  }
  
  all_results[[name]] <- plot_nmda_scatter(
    rdata_file = a$rdata_file,
    output_prefix = a$output_prefix,
    plot_title = a$plot_title,
    x_label = a$x_label,
    mode_idx = 2,
    flip_sign = TRUE
  )
  
  print(all_results[[name]]$plot)
}

# ------------------------------------------------------------------
# 5. Save combined r summary
# ------------------------------------------------------------------

r_summary <- do.call(
  rbind,
  lapply(names(all_results), function(name) {
    all_results[[name]]$summary
  })
)

write.csv(
  r_summary,
  file.path(out_dir, "NMDA_scatter_r_summary.csv"),
  row.names = FALSE
)

print(r_summary)

cat("\nNMDA scatter plots and source data saved.\n")