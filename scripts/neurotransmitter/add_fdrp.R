# -------------------------------------------------------------------------
# Compute FDR-corrected p-values for neurotransmitter regression results
# across all datasets: CMI-math, CMI-reading, Stanford-math, Stanford-reading
# Author: Yuan Zhang
# Date: 2025-07-25
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
  df$fdrp <- 1
  
  # Define age bins present in brain_map column
  agebins <- c("Full", "7-8.9y", "9-10.9y", "11-13.9y")
  
  # Compute FDR correction for each age bin separately
  for (k in 1:length(agebins)) {
    idx <- which(grepl(agebins[k], df$brain_map))
    df$fdrp[idx] <- p.adjust(df$moran_p[idx], method = "fdr", n = length(idx))
  }
  
  # Save the updated DataFrame with FDR-corrected p-values
  write.csv(df, output_file, row.names = FALSE)
  
  cat("FDR-corrected results saved to:", output_file, "\n")
}

# -------------------------------------------------------------------------
# Define datasets to process
# -------------------------------------------------------------------------
datasets <- list(
  list(
    input = "results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_math.csv",
    output = "results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_math_wFDRp.csv"
  ),
  list(
    input = "results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_read.csv",
    output = "results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_read_wFDRp.csv"
  ),
  list(
    input = "results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_math.csv",
    output = "results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_math_wFDRp.csv"
  ),
  list(
    input = "results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_read.csv",
    output = "results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_read_wFDRp.csv"
  )
)

# -------------------------------------------------------------------------
# Run FDR computation for all datasets
# -------------------------------------------------------------------------
for (ds in datasets) {
  compute_fdr(ds$input, ds$output)
}
