# -------------------------------------------------------------------------
# Author: Yuan Zhang
# Date: 2025-07-25
# Compute Bayes Factors (BF10) for association between brain GMV weight maps 
# (from CCA) and neurotransmitter receptor maps for CMI datasets:
#   - CMI Math
#   - CMI Reading
# -------------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(BayesFactor)

# -------------------------------------------------------------------------
# Load neurotransmitter receptor data
# -------------------------------------------------------------------------
fname_receptors <- 'data/neurotransmitter/receptor_data_bn246.csv'
receptors <- read.csv(fname_receptors, header = FALSE)

# Receptor names (correspond to columns in receptors CSV)
receptor_names <- c(
  '5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
  'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
  'NMDA', 'VAChT'
)

# -------------------------------------------------------------------------
# Function to compute Bayes Factors for a given brain map dataset
# -------------------------------------------------------------------------
compute_bayes_factors <- function(bmap_file, bmap_name, output_prefix) {
  # Load GMV weight maps for all age bins
  bmaps <- read.csv(bmap_file)
  
  # Combine GMV maps with receptor data (only first 218 rows are used)
  df <- data.frame(bmaps, receptors[1:218, ])
  
  # Assign receptor column names (starting from column 5 onward)
  colnames(df)[5:ncol(df)] <- receptor_names
  
  # Ensure valid column names for formulas
  colnames(df) <- make.names(colnames(df))
  
  # Age bins (columns 1–4 in bmaps)
  agebins <- c("Full", "7-8.9y", "9-10.9y", "11-13.9y")
  
  # Loop through each age bin
  for (k in 1:length(agebins)) {
    # Output file for this age bin
    outputf <- paste(output_prefix, 'bf10_', bmap_name, '_receptors_', agebins[k], '.csv', sep = "")
    
    bf10 <- c()
    
    # Loop through each receptor
    for (i in 5:ncol(df)) {
      # Construct formula: brain map ~ receptor
      formula <- as.formula(paste(colnames(df)[k], " ~ ", colnames(df)[i], sep = ""))
      
      # Compute Bayes Factor (BF10)
      bf <- lmBF(formula, data = df)
      bf_value <- exp(bf@bayesFactor$bf)
      bf10 <- c(bf10, bf_value)
    }
    
    # Create results dataframe
    res <- data.frame(receptor = receptor_names, BF10 = bf10)
    res$bmap <- bmap_name
    res$agbin <- agebins[k]
    
    # Reorder rows (specific receptor order)
    idx <- c(9:11, 15, 18, 12, 1:7, 14, 19, 17, 8, 13, 16)
    res <- res[idx, c(3, 1, 2, 4)]
    
    # Save results
    write.csv(res, file = outputf, row.names = FALSE)
    cat("Bayes Factor results saved to:", outputf, "\n")
  }
}

# -------------------------------------------------------------------------
# Run for both CMI Math and CMI Reading
# -------------------------------------------------------------------------
# CMI Math
compute_bayes_factors(
  bmap_file = 'results/age_analysis/cmi/gmv_weights_agebin_math.csv',
  bmap_name = 'math_gmv',
  output_prefix = 'results/neurotransmitter/cmi/'
)

# CMI Reading
compute_bayes_factors(
  bmap_file = 'results/age_analysis/cmi/gmv_weights_agebin_read.csv',
  bmap_name = 'read_gmv',
  output_prefix = 'results/neurotransmitter/cmi/'
)
