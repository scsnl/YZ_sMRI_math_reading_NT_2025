# -------------------------------------------------------------------------
# Scatter Plot Analysis of GMV Weights vs. Neurotransmitter Receptor Maps
# Author: Yuan Zhang
# Date: 2025-07-25
# 
# This script:
#   - Loads GMV weights (from CCA) for both CMI and Stanford datasets.
#   - Plots GMV vs. receptor densities for both math and reading datasets.
#   - Generates scatter plots for a specified receptor (or all receptors).
#
# Datasets processed:
#   - CMI math
#   - Stanford math
#   - CMI reading
#   - Stanford reading
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ggplot2)

# -------------------------------------------------------------------------
# Load Neurotransmitter Receptor Data
# -------------------------------------------------------------------------
fname_receptors <- 'data/neurotransmitter/receptor_data_bn246.csv'
receptors <- read.csv(fname_receptors, header = FALSE)

# Names of the 19 receptor types
receptor_names <- c(
  '5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
  'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
  'NMDA', 'VAChT'
)

# -------------------------------------------------------------------------
# Function to generate scatter plots for one dataset (CMI or Stanford)
# -------------------------------------------------------------------------
plot_scatter <- function(bmap_file, dataset_name, bmap_name, receptor_idx) {
  # Load GMV weight data
  bmaps <- read.csv(bmap_file)
  
  # Combine GMV map (Full sample) with receptor data
  df <- data.frame(bmaps$Full, receptors[1:218, ])
  colnames(df) <- c('bmap', receptor_names)
  
  # Get receptor name
  receptor <- receptor_names[receptor_idx]
  
  # Output file path
  outputf <- paste0('results/network_analysis/', dataset_name, '_scatter_', bmap_name, '_', receptor, '.eps')
  
  # Plot GMV vs. receptor density
  ggplot(df, aes(x = bmap, y = df[, receptor_idx + 1])) + 
    geom_point(size = 3, alpha = 0.5) +   
    ylab(receptor) +
    geom_smooth(method = lm, se = TRUE, color = "black") +
    theme_classic() +
    theme(legend.position = "none")
  
  # Save plot
  ggsave(outputf, units = "in", width = 5, height = 5, device = cairo_ps)
  
  cat("Scatter plot saved to:", outputf, "\n")
}

# -------------------------------------------------------------------------
# Generate Plots for Math and Reading (CMI and Stanford)
# -------------------------------------------------------------------------
i <- 18  # Example: NMDA receptor (index 18)

# Math datasets
plot_scatter('results/age_analysis/cmi/gmv_weights_agebin_math.csv', 'cmi', 'math_gmv', i)
plot_scatter('results/age_analysis/stanford/gmv_weights_agebin_math.csv', 'stanford', 'math_gmv', i)

# Reading datasets
plot_scatter('results/age_analysis/cmi/gmv_weights_agebin_read.csv', 'cmi', 'read_gmv', i)
plot_scatter('results/age_analysis/stanford/gmv_weights_agebin_read.csv', 'stanford', 'read_gmv', i)

# -------------------------------------------------------------------------
# To generate plots for ALL receptors, you can do:
# for (i in 1:19) {
#   plot_scatter('results/age_analysis/cmi/gmv_weights_agebin_math.csv', 'cmi', 'math_gmv', i)
#   plot_scatter('results/age_analysis/stanford/gmv_weights_agebin_math.csv', 'stanford', 'math_gmv', i)
#   plot_scatter('results/age_analysis/cmi/gmv_weights_agebin_read.csv', 'cmi', 'read_gmv', i)
#   plot_scatter('results/age_analysis/stanford/gmv_weights_agebin_read.csv', 'stanford', 'read_gmv', i)
# }
# -------------------------------------------------------------------------
