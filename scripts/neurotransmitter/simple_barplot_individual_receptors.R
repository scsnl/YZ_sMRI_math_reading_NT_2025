############################################################
# Visualization of Adjusted R² for Neurotransmitter Analysis
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script generates horizontal barplots of adjusted R² values 
# for the relationship between GMV CCA modes and neurotransmitter 
# receptor maps across different age bins.
#
# It is configured to handle:
#   - CMI math
#   - CMI reading
#   - Stanford math
#   - Stanford reading
#
# The output plots are saved as EPS files for publication-quality figures.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory to project root
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear the environment
rm(list = ls())

# Load required libraries
library(ggplot2)    # For plotting
library(forcats)    # For factor reordering

# ------------------------------------------------------------------
# 2. Define Cohorts and Tasks
# ------------------------------------------------------------------

# Each entry in the list specifies:
#  - Path to the input CSV (individual regression results)
#  - Brain map name (used in plot titles/filenames)
#  - Cohort name (CMI or Stanford)
tasks <- list(
  cmi_reading = list(
    fname   = 'results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_read.csv',
    bmap    = "read_gmv",
    cohort  = "cmi",
    outdir  = "results/neurotransmitter/cmi/"
  ),
  cmi_math = list(
    fname   = 'results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_math.csv',
    bmap    = "math_gmv",
    cohort  = "cmi",
    outdir  = "results/neurotransmitter/cmi/"
  ),
  stanford_reading = list(
    fname   = 'results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_read.csv',
    bmap    = "read_gmv",
    cohort  = "stanford",
    outdir  = "results/neurotransmitter/stanford/"
  ),
  stanford_math = list(
    fname   = 'results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_math.csv',
    bmap    = "math_gmv",
    cohort  = "stanford",
    outdir  = "results/neurotransmitter/stanford/"
  )
)

# Define age bins
agebins <- c("Full", "7-8.9y", "9-10.9y", "11-13.9y")

# ------------------------------------------------------------------
# 3. Generate Barplots for Each Task
# ------------------------------------------------------------------

for (task_name in names(tasks)) {
  cat("\nProcessing task:", task_name, "\n")
  
  # Extract task parameters
  fname     <- tasks[[task_name]]$fname
  bmap      <- tasks[[task_name]]$bmap
  cohort    <- tasks[[task_name]]$cohort
  output_path <- tasks[[task_name]]$outdir
  
  # Read the regression results CSV
  df_all <- read.csv(fname)
  
  # Loop through each age bin
  for (k in 1:length(agebins)) {
    agebin <- agebins[k]
    cat("  - Generating barplot for age bin:", agebin, "\n")
    
    # Filter rows corresponding to the current age bin
    df <- df_all[grepl(agebin, df_all$brain_map), ]
    
    # Define output filename
    postfix <- paste("barplot_adj_r2_", bmap, "_individual_neurotransmitter_",
                     cohort, "_", agebin, sep = "")
    outputf <- paste(output_path, '/', postfix, '.eps', sep = "")
    
    # Set receptor_set as factor with custom order for consistent plots
    df$receptor_set <- as.factor(df$receptor_set)
    df$receptor_set <- factor(df$receptor_set,
                              levels = levels(df$receptor_set)[c(9:11, 15, 18, 12, 1:6, 7, 14, 19, 17, 8, 13, 16)])
    
    # Generate horizontal barplot
    p <- ggplot(df, aes(x = fct_rev(receptor_set), y = adj_r2)) + 
      geom_bar(stat = "identity", width = 0.86, fill = "darkorange") +
      ylim(c(-0.01, 0.6)) +
      coord_flip() +
      theme_classic() +
      theme(
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")
      )
    
    # Save the plot
    ggsave(outputf, plot = p, units = "in", width = 5, height = 4)
  }
}

cat("\nAll visualizations completed!\n")
