# -------------------------------------------------------------------------
# Network Analysis: Correlation of GMV maps with NMDA receptor density within
# each of Shirer networks
# Author: Yuan Zhang
# Date: 2025-07-25
# 
# This script:
#   - Loads GMV weight maps (from CCA) for both CMI and Stanford cohorts.
#   - Runs correlation analysis between GMV maps and NMDA receptor densities
#     for each Shirer functional network.
#   - Saves results (t-statistic, r, p-value, FDR-corrected p) and plots.
#
# Cohorts processed:
#   - CMI (math, reading)
#   - Stanford (math, reading)
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ggplot2)

# -------------------------------------------------------------------------
# Load Shirer network mapping
# -------------------------------------------------------------------------
bn_shirer <- read.csv("data/atlas/shirer_to_bn.csv", header = FALSE)
colnames(bn_shirer) <- c("BN_ID", "Shirer_ID", "Shirer_network")
bn_shirer <- bn_shirer[1:218, ]  # Keep only 218 cortical ROIs
bn_shirer$Shirer_network <- trimws(bn_shirer$Shirer_network, which = "left")

# -------------------------------------------------------------------------
# Load neurotransmitter receptor data
# -------------------------------------------------------------------------
fname_receptors <- "data/neurotransmitter/receptor_data_bn246.csv"
receptors <- read.csv(fname_receptors, header = FALSE)

receptor_names <- c(
  '5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
  'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
  'NMDA', 'VAChT'
)

# -------------------------------------------------------------------------
# Function to run network-level correlation analysis for NMDA
# -------------------------------------------------------------------------
analyze_network <- function(bmap_file, cohort_name, bmap_name) {
  # Load GMV weight maps
  bmaps <- read.csv(bmap_file)
  
  # Extract the first column (Full sample GMV map)
  bmap_full <- bmaps$Full
  
  # Combine GMV weights (Full sample) and NMDA receptor density
  df <- data.frame(bmap = bmap_full, NMDA = receptors[1:218, which(receptor_names == "NMDA")])
  df$Shirer_network <- bn_shirer$Shirer_network
  
  # Prepare result data frame
  networks <- unique(df$Shirer_network)
  nn <- length(networks)
  res <- data.frame(network = networks, t = numeric(nn), df = numeric(nn),
                    r = numeric(nn), p = numeric(nn))
  
  # Create output folder if not exists
  output_dir <- paste0("results/network_analysis/", cohort_name, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Correlation analysis for each Shirer network
  for (j in 1:nn) {
    df_net <- subset(df, Shirer_network == networks[j])
    tmp <- cor.test(df_net$bmap, df_net$NMDA)
    
    res$t[j] <- tmp$statistic
    res$df[j] <- tmp$parameter
    res$r[j] <- tmp$estimate
    res$p[j] <- tmp$p.value
    
    # Plot GMV vs NMDA for this network
    ylabel <- paste("NMDA (", networks[j], ")", sep = "")
    xlabel <- paste("Brain map of ", bmap_name, " (", networks[j], ")", sep = "")
    
    p_fig <- ggplot(df_net, aes(x = bmap, y = NMDA)) +
      geom_point(size = 5, alpha = 0.5) +
      ylab(ylabel) +
      xlab(xlabel) +
      geom_smooth(method = lm, se = TRUE) +
      theme_classic() +
      theme(legend.position = "none")
    
    # Save plot
    output_fig <- paste0(output_dir, "nmda_", bmap_name, "_", networks[j], ".eps")
    ggsave(output_fig, p_fig, units = "in", width = 5, height = 5, device = cairo_ps)
  }
  
  # Add FDR-corrected p-values
  res$fdrp <- p.adjust(res$p, method = "fdr", n = length(res$p))
  
  # Save CSV results
  output_csv <- paste0(output_dir, "nmda_", bmap_name, "_shirer.csv")
  write.csv(res, file = output_csv, row.names = FALSE)
  
  cat("Analysis completed for:", cohort_name, "-", bmap_name, "\n")
}

# -------------------------------------------------------------------------
# Run for both math and reading for CMI and Stanford
# -------------------------------------------------------------------------
# CMI
analyze_network("results/age_analysis/cmi/gmv_weights_agebin_math.csv", "cmi", "math_gmv")
analyze_network("results/age_analysis/cmi/gmv_weights_agebin_read.csv", "cmi", "read_gmv")

# Stanford
analyze_network("results/age_analysis/stanford/gmv_weights_agebin_math.csv", "stanford", "math_gmv")
analyze_network("results/age_analysis/stanford/gmv_weights_agebin_read.csv", "stanford", "read_gmv")