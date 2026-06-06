# -------------------------------------------------------------------------
# Network Analysis: Correlation of GMV maps with NMDA receptor density within
# each Shirer network
# Author: Yuan Zhang
# Date: 2026-05-29
#
# This script:
#   - Loads CCA .RData files for CMI and Stanford cohorts.
#   - Extracts sign-oriented GMV structural phenotype maps from CCA mode 2.
#   - Runs Pearson correlation analyses between GMV maps and NMDA receptor
#     density within each Shirer functional network.
#   - Saves t statistic, df, r, p value, 95% CI, FDR-corrected p value,
#     full reporting strings, plots, and source data.
#
# Cohorts processed:
#   - CMI math
#   - CMI reading
#   - Stanford math
#   - Stanford reading
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ggplot2)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

format_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  } else if (p < 0.001) {
    return(formatC(p, format = "e", digits = 2))
  } else {
    return(sprintf("%.3f", p))
  }
}

format_ci <- function(low, high) {
  sprintf("[%.2f, %.2f]", low, high)
}

# Extract sign-oriented GMV map from CCA .RData file.
# This computes correlations between original ROI-level GMV values and
# the brain canonical variate for the selected mode.
get_mode_gmv_weights <- function(rdata_file, mode_idx = 2, flip_sign = TRUE) {
  load(rdata_file)
  
  bmap <- as.vector(cor(Xorig, res$Cx[, mode_idx]))
  
  if (flip_sign) {
    bmap <- -1 * bmap
  }
  
  return(bmap)
}

# Safely run Pearson correlation.
safe_cor_test <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])
  
  if (length(x) < 3) {
    return(NULL)
  }
  
  if (sd(x) == 0 || sd(y) == 0) {
    return(NULL)
  }
  
  cor.test(x, y, method = "pearson")
}

# -------------------------------------------------------------------------
# Load Shirer network mapping
# -------------------------------------------------------------------------

bn_shirer <- read.csv("data/atlas/shirer_to_bn.csv", header = FALSE)
colnames(bn_shirer) <- c("BN_ID", "Shirer_ID", "Shirer_network")
bn_shirer <- bn_shirer[1:218, ]
bn_shirer$Shirer_network <- trimws(bn_shirer$Shirer_network)

# -------------------------------------------------------------------------
# Load neurotransmitter receptor data
# -------------------------------------------------------------------------

fname_receptors <- "data/neurotransmitter/receptor_data_bn246.csv"
receptors <- read.csv(fname_receptors, header = FALSE)

receptor_names <- c(
  "5HT1a", "5HT1b", "5HT2a", "5HT4", "5HT6", "5HTT", "A4B2", "CB1",
  "D1", "D2", "DAT", "GABAa", "H3", "M1", "mGluR5", "MOR", "NET",
  "NMDA", "VAChT"
)

nmda_map <- as.numeric(receptors[1:218, which(receptor_names == "NMDA")])

# -------------------------------------------------------------------------
# Function to run network-level correlation analysis for NMDA
# -------------------------------------------------------------------------

analyze_network <- function(
    rdata_file,
    cohort_name,
    bmap_name,
    mode_idx = 2,
    flip_sign = TRUE
) {
  
  cat("\nRunning network analysis for:", cohort_name, "-", bmap_name, "\n")
  
  # -----------------------------------------------------------------------
  # Extract sign-oriented GMV structural phenotype from CCA .RData
  # -----------------------------------------------------------------------
  
  bmap_full <- get_mode_gmv_weights(
    rdata_file = rdata_file,
    mode_idx = mode_idx,
    flip_sign = flip_sign
  )
  
  # Combine GMV weights, NMDA receptor density, and network labels
  df <- data.frame(
    ROI = seq_len(218),
    BN_ID = bn_shirer$BN_ID,
    Shirer_network = bn_shirer$Shirer_network,
    bmap = as.numeric(bmap_full),
    NMDA = nmda_map
  )
  
  df <- df[complete.cases(df), ]
  
  # Output folder
  output_dir <- paste0("results/network_analysis/", cohort_name, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save full ROI-level source data used for all network analyses
  source_data_file <- paste0(
    output_dir,
    "source_data_nmda_",
    bmap_name,
    "_all_networks.csv"
  )
  
  write.csv(df, source_data_file, row.names = FALSE)
  
  # Networks
  networks <- unique(df$Shirer_network)
  networks <- networks[!is.na(networks)]
  
  # Results data frame
  res <- data.frame(
    network = character(),
    n_ROI = numeric(),
    t = numeric(),
    df = numeric(),
    r = numeric(),
    CI_low = numeric(),
    CI_high = numeric(),
    p = numeric(),
    stringsAsFactors = FALSE
  )
  
  # -----------------------------------------------------------------------
  # Correlation analysis for each Shirer network
  # -----------------------------------------------------------------------
  
  for (j in seq_along(networks)) {
    
    net_name <- networks[j]
    df_net <- subset(df, Shirer_network == net_name)
    
    tmp <- safe_cor_test(df_net$bmap, df_net$NMDA)
    
    if (is.null(tmp)) {
      tmp_res <- data.frame(
        network = net_name,
        n_ROI = nrow(df_net),
        t = NA,
        df = NA,
        r = NA,
        CI_low = NA,
        CI_high = NA,
        p = NA,
        stringsAsFactors = FALSE
      )
    } else {
      tmp_res <- data.frame(
        network = net_name,
        n_ROI = nrow(df_net),
        t = unname(tmp$statistic),
        df = unname(tmp$parameter),
        r = unname(tmp$estimate),
        CI_low = tmp$conf.int[1],
        CI_high = tmp$conf.int[2],
        p = tmp$p.value,
        stringsAsFactors = FALSE
      )
    }
    
    res <- rbind(res, tmp_res)
    
    # ---------------------------------------------------------------------
    # Save network-specific source data
    # ---------------------------------------------------------------------
    
    safe_net_name <- gsub("[^A-Za-z0-9_]+", "_", net_name)
    
    write.csv(
      df_net,
      paste0(
        output_dir,
        "source_data_nmda_",
        bmap_name,
        "_",
        safe_net_name,
        ".csv"
      ),
      row.names = FALSE
    )
    
    # ---------------------------------------------------------------------
    # Plot GMV vs NMDA for this network
    # ---------------------------------------------------------------------
    
    ylabel <- paste0("NMDA (", net_name, ")")
    xlabel <- paste0("Brain map of ", bmap_name, " (", net_name, ")")
    
    p_fig <- ggplot(df_net, aes(x = bmap, y = NMDA)) +
      geom_point(size = 5, alpha = 0.5) +
      geom_smooth(method = lm, se = TRUE) +
      ylab(ylabel) +
      xlab(xlabel) +
      theme_classic(base_size = 14) +
      theme(legend.position = "none")
    
    output_fig_eps <- paste0(
      output_dir,
      "nmda_",
      bmap_name,
      "_",
      safe_net_name,
      ".eps"
    )
    
    output_fig_pdf <- paste0(
      output_dir,
      "nmda_",
      bmap_name,
      "_",
      safe_net_name,
      ".pdf"
    )
    
    output_fig_png <- paste0(
      output_dir,
      "nmda_",
      bmap_name,
      "_",
      safe_net_name,
      ".png"
    )
    
    ggsave(output_fig_eps, p_fig, units = "in", width = 5, height = 5, device = cairo_ps)
    ggsave(output_fig_pdf, p_fig, units = "in", width = 5, height = 5)
    ggsave(output_fig_png, p_fig, units = "in", width = 5, height = 5, dpi = 300)
    
    # Save fitted regression line and 95% CI band source data
    if (nrow(df_net) >= 3 && sd(df_net$bmap, na.rm = TRUE) > 0 && sd(df_net$NMDA, na.rm = TRUE) > 0) {
      
      lm_fit <- lm(NMDA ~ bmap, data = df_net)
      
      x_grid <- data.frame(
        bmap = seq(
          min(df_net$bmap, na.rm = TRUE),
          max(df_net$bmap, na.rm = TRUE),
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
        bmap = x_grid$bmap,
        fitted_NMDA = pred_ci[, "fit"],
        ci_lower = pred_ci[, "lwr"],
        ci_upper = pred_ci[, "upr"]
      )
      
      write.csv(
        df_fit_ci,
        paste0(
          output_dir,
          "source_data_nmda_",
          bmap_name,
          "_",
          safe_net_name,
          "_fit_95CI.csv"
        ),
        row.names = FALSE
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # FDR correction across networks
  # -----------------------------------------------------------------------
  
  res$fdrp <- NA
  valid_p <- !is.na(res$p)
  res$fdrp[valid_p] <- p.adjust(res$p[valid_p], method = "fdr")
  
  # -----------------------------------------------------------------------
  # Reporting strings
  # -----------------------------------------------------------------------
  
  res$Full_report <- mapply(
    function(t, df, r, ci_low, ci_high, p) {
      if (any(is.na(c(t, df, r, ci_low, ci_high, p)))) {
        return(NA_character_)
      }
      
      sprintf(
        "t(%d) = %.2f, p = %s, r = %.2f, 95%% CI %s",
        round(df),
        t,
        format_p(p),
        r,
        format_ci(ci_low, ci_high)
      )
    },
    res$t,
    res$df,
    res$r,
    res$CI_low,
    res$CI_high,
    res$p
  )
  
  res$Full_report_FDR <- mapply(
    function(report, fdrp) {
      if (is.na(report) || is.na(fdrp)) {
        return(NA_character_)
      }
      
      sprintf(
        "%s, FDR-corrected p = %s",
        report,
        format_p(fdrp)
      )
    },
    res$Full_report,
    res$fdrp
  )
  
  # -----------------------------------------------------------------------
  # Save CSV results
  # -----------------------------------------------------------------------
  
  output_csv <- paste0(output_dir, "nmda_", bmap_name, "_shirer.csv")
  write.csv(res, file = output_csv, row.names = FALSE)
  
  cat("Analysis completed for:", cohort_name, "-", bmap_name, "\n")
  cat("Results saved to:", output_csv, "\n")
  
  return(res)
}

# -------------------------------------------------------------------------
# Run for CMI and Stanford, math and reading
# -------------------------------------------------------------------------

results_all <- list()

results_all$cmi_math <- analyze_network(
  rdata_file = "results/cca/cmi/wholebrain_cca_cmi_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
  cohort_name = "cmi",
  bmap_name = "math_gmv",
  mode_idx = 2,
  flip_sign = TRUE
)

results_all$cmi_reading <- analyze_network(
  rdata_file = "results/cca/cmi/wholebrain_cca_cmi_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
  cohort_name = "cmi",
  bmap_name = "read_gmv",
  mode_idx = 2,
  flip_sign = TRUE
)

results_all$stanford_math <- analyze_network(
  rdata_file = "results/cca/stanford/wholebrain_cca_stanford_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
  cohort_name = "stanford",
  bmap_name = "math_gmv",
  mode_idx = 2,
  flip_sign = TRUE
)

results_all$stanford_reading <- analyze_network(
  rdata_file = "results/cca/stanford/wholebrain_cca_stanford_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
  cohort_name = "stanford",
  bmap_name = "read_gmv",
  mode_idx = 2,
  flip_sign = TRUE
)

# -------------------------------------------------------------------------
# Save combined result table
# -------------------------------------------------------------------------

combined_results <- do.call(
  rbind,
  lapply(names(results_all), function(name) {
    tmp <- results_all[[name]]
    tmp$analysis <- name
    tmp
  })
)

combined_results <- combined_results[, c(
  "analysis",
  "network",
  "n_ROI",
  "t",
  "df",
  "r",
  "CI_low",
  "CI_high",
  "p",
  "fdrp",
  "Full_report",
  "Full_report_FDR"
)]

dir.create("results/network_analysis/combined", showWarnings = FALSE, recursive = TRUE)

write.csv(
  combined_results,
  "results/network_analysis/combined/nmda_network_analysis_all_results.csv",
  row.names = FALSE
)

print(combined_results)