# -------------------------------------------------------------------------
# Compute FDR-corrected p-values for neurotransmitter regression results
# for joint CCA Mode 2 analyses
# Author: Yuan Zhang
# Date: 2026-04-13
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

# -------------------------------------------------------------------------
# Function to compute FDR p-values and update CSV
# -------------------------------------------------------------------------
compute_fdr <- function(input_file, output_file) {
  # Load CSV with original p-values
  df <- read.csv(input_file)
  
  # Initialize a new column for FDR-corrected p-values
  df$fdrp <- p.adjust(df$moran_p, method = "fdr")
  
  # Save the updated DataFrame with FDR-corrected p-values
  write.csv(df, output_file, row.names = FALSE)
  
  cat("FDR-corrected results saved to:", output_file, "\n")
}

# -------------------------------------------------------------------------
# Define datasets to process
# -------------------------------------------------------------------------
datasets <- list(
  list(
    input  = "results/neurotransmitter/cmi/shared/individual_neurotransmitter_regression_results_shared_mode2.csv",
    output = "results/neurotransmitter/cmi/shared/individual_neurotransmitter_regression_results_shared_mode2_wFDRp.csv"
  ),
  list(
    input  = "results/neurotransmitter/stanford/shared/individual_neurotransmitter_regression_results_shared_mode2.csv",
    output = "results/neurotransmitter/stanford/shared/individual_neurotransmitter_regression_results_shared_mode2_wFDRp.csv"
  )
)

# -------------------------------------------------------------------------
# Run FDR computation for all datasets
# -------------------------------------------------------------------------
for (ds in datasets) {
  compute_fdr(ds$input, ds$output)
}