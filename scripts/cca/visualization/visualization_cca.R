############################################################
# Visualize res$corr.Y.Cy as heatmaps
# for CMI and Stanford combined CCA analyses
# (flip all signs before plotting)
#
# Also create a scatter plot for a specific mode
# (default: Mode 2) using:
#   brain = res$Cx[, mode]
#   beh   = res$Cy[, mode]
#
# Saves:
#   - CSV of flipped matrices
#   - heatmap: PNG / PDF / EPS
#   - scatter: EPS (cairo_ps) / PNG / PDF
#
# Author: Yuan Zhang
# Date: 2026-04-03
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ggplot2)
library(reshape2)

# ------------------------------------------------------------------
# 2. Helper function
# ------------------------------------------------------------------
plot_corr_y_cy_heatmap_and_scatter <- function(
    rdata_file,
    output_dir,
    cohort_label,
    mode_to_plot = 2
) {
  
  # Load result
  load(rdata_file)
  
  # -----------------------------
  # A. Heatmap from corr.Y.Cy
  # -----------------------------
  ycy_mat <- res$corr.Y.Cy
  ycy_mat_flip <- -1 * ycy_mat
  
  pretty_rows <- c(
    "Age",
    "Numerical Operation",
    "Math Problem Solving",
    "Word Reading",
    "Reading Comprehension"
  )
  
  if (nrow(ycy_mat_flip) == length(pretty_rows)) {
    rownames(ycy_mat_flip) <- pretty_rows
  }
  
  colnames(ycy_mat_flip) <- paste0("Mode", seq_len(ncol(ycy_mat_flip)))
  
  write.csv(
    ycy_mat_flip,
    file.path(output_dir, "corr_Y_Cy_flipped.csv"),
    row.names = TRUE
  )
  
  plot_df <- reshape2::melt(ycy_mat_flip)
  colnames(plot_df) <- c("Behavior", "Mode", "Correlation")
  
  plot_df$Behavior <- factor(
    plot_df$Behavior,
    levels = rev(rownames(ycy_mat_flip))
  )
  
  plot_df$Mode <- factor(
    plot_df$Mode,
    levels = colnames(ycy_mat_flip)
  )
  
  p_heat <- ggplot(plot_df, aes(x = Mode, y = Behavior, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 5) +
    scale_fill_gradient2(
      low = "steelblue",
      mid = "white",
      high = "orange",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    labs(
      title = paste0(cohort_label, ": Correlations Between Behavioral Variables and Canonical Variates"),
      x = "",
      y = "",
      fill = "r"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave(
    file.path(output_dir, "corr_Y_Cy_flipped_heatmap.png"),
    p_heat,
    width = 8,
    height = 4.8,
    dpi = 300
  )
  
  ggsave(
    file.path(output_dir, "corr_Y_Cy_flipped_heatmap.pdf"),
    p_heat,
    width = 8,
    height = 4.8
  )
  
  postscript(
    file = file.path(output_dir, "corr_Y_Cy_flipped_heatmap.eps"),
    width = 8,
    height = 4.8,
    horizontal = FALSE,
    onefile = FALSE,
    paper = "special"
  )
  print(p_heat)
  dev.off()
  
  # -----------------------------
  # B. Scatter plot for one mode
  # -----------------------------
  stopifnot(mode_to_plot >= 1, mode_to_plot <= ncol(res$Cx), mode_to_plot <= ncol(res$Cy))
  
  df_scatter <- data.frame(
    brain = res$Cx[, mode_to_plot],
    beh   = res$Cy[, mode_to_plot]
  )
  
  p_scatter <- ggplot(df_scatter, aes(x = brain, y = beh)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = lm, se = TRUE) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0(cohort_label, ": Brain vs Behavioral Canonical Variates (Mode ", mode_to_plot, ")"),
      x = paste0("Brain canonical variate (Mode ", mode_to_plot, ")"),
      y = paste0("Behavioral canonical variate (Mode ", mode_to_plot, ")")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  # EPS editable in Illustrator
  ggsave(
    filename = file.path(output_dir, paste0("scatter_mode", mode_to_plot, ".eps")),
    plot = p_scatter,
    units = "in",
    width = 5,
    height = 5,
    device = cairo_ps
  )
  
  # Optional additional versions
  ggsave(
    filename = file.path(output_dir, paste0("scatter_mode", mode_to_plot, ".png")),
    plot = p_scatter,
    units = "in",
    width = 5,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0("scatter_mode", mode_to_plot, ".pdf")),
    plot = p_scatter,
    units = "in",
    width = 5,
    height = 5
  )
  
  return(list(heatmap = p_heat, scatter = p_scatter))
}

# ------------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------------

# CMI
cmi_rdata <- paste0(
  "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/",
  "CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)
cmi_outdir <- "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/"

# Stanford
stan_rdata <- paste0(
  "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined/",
  "CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
)
stan_outdir <- "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined/"

# ------------------------------------------------------------------
# 4. Run for both cohorts
# ------------------------------------------------------------------
plots_cmi <- plot_corr_y_cy_heatmap_and_scatter(
  rdata_file = cmi_rdata,
  output_dir = cmi_outdir,
  cohort_label = "CMI Combined CCA",
  mode_to_plot = 2
)

plots_stan <- plot_corr_y_cy_heatmap_and_scatter(
  rdata_file = stan_rdata,
  output_dir = stan_outdir,
  cohort_label = "Stanford Combined CCA",
  mode_to_plot = 2
)

print(plots_cmi$heatmap)
print(plots_cmi$scatter)
print(plots_stan$heatmap)
print(plots_stan$scatter)