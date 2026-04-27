# -------------------------------------------------------------------------
# Author: Yuan Zhang
# Date: 2026-04-13
#
# Compute Bayes Factors (BF10) for association between Mode 2 brain GMV
# weight maps and neurotransmitter receptor maps for CMI joint CCA model  
#
# Brain map source:
#   *_coef.csv
# Weight column used:
#   xloading_m2
# Sign is flipped by multiplying by -1 before analysis
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(BayesFactor)

# -------------------------------------------------------------------------
# Load neurotransmitter receptor data
# -------------------------------------------------------------------------
fname_receptors <- "data/neurotransmitter/receptor_data_bn246.csv"
receptors <- read.csv(fname_receptors, header = FALSE)

# Receptor names (correspond to columns in receptors CSV)
receptor_names <- c(
  "5HT1a", "5HT1b", "5HT2a", "5HT4", "5HT6", "5HTT", "A4B2", "CB1",
  "D1", "D2", "DAT", "GABAa", "H3", "M1", "mGluR5", "MOR", "NET",
  "NMDA", "VAChT"
)

# -------------------------------------------------------------------------
# Function to compute Bayes Factors for a single Mode 2 brain map
# -------------------------------------------------------------------------
compute_bayes_factors_mode2 <- function(coef_file, bmap_name, output_file,
                                        weight_col = "xloading_m2",
                                        flip_sign = TRUE) {
  # Load coefficient file
  coef_df <- read.csv(coef_file)
  
  if (!(weight_col %in% colnames(coef_df))) {
    stop(paste(weight_col, "not found in", coef_file))
  }
  
  # Extract Mode 2 brain map
  bmap <- coef_df[[weight_col]]
  
  # Flip sign if requested
  if (flip_sign) {
    bmap <- -1 * bmap
  }
  
  # Use only first 218 rows to match receptor maps
  bmap <- bmap[1:218]
  
  # Build dataframe
  df <- data.frame(
    bmap = bmap,
    receptors[1:218, ]
  )
  
  # Assign receptor names
  colnames(df)[2:ncol(df)] <- receptor_names
  
  # Ensure valid column names for formulas
  colnames(df) <- make.names(colnames(df))
  
  # Compute BF10 for each receptor
  bf10 <- c()
  
  for (i in 2:ncol(df)) {
    formula <- as.formula(paste("bmap ~", colnames(df)[i]))
    bf <- lmBF(formula, data = df)
    bf_value <- exp(bf@bayesFactor$bf)
    bf10 <- c(bf10, bf_value)
  }
  
  # Create results dataframe
  res <- data.frame(
    bmap = bmap_name,
    receptor = receptor_names,
    BF10 = bf10
  )
  
  # Reorder rows to match previous output style
  idx <- c(9:11, 15, 18, 12, 1:7, 14, 19, 17, 8, 13, 16)
  res <- res[idx, c("bmap", "receptor", "BF10")]
  
  # Save results
  write.csv(res, file = output_file, row.names = FALSE)
  cat("Bayes Factor results saved to:", output_file, "\n")
}

# -------------------------------------------------------------------------
# Run for current CMI domain-specific datasets
# -------------------------------------------------------------------------

# CMI shared
compute_bayes_factors_mode2(
  coef_file = "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_coef.csv",
  bmap_name = "CMI_shared_mode2",
  output_file = "results/neurotransmitter/cmi/shared/bf10_CMI_shared_mode2_receptors.csv"
)
