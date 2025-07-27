############################################################
# Prediction Using CCA Models (CMI → Stanford)
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script applies Canonical Correlation Analysis (CCA) models
# derived from the CMI cohort to predict brain–behavior scores 
# in the Stanford cohort (both math and reading models).
#
# Steps:
#  1. Load Stanford behavioral and brain data.
#  2. Load the CCA model trained on CMI data (math and reading).
#  3. Project Stanford brain data onto the PCA space derived from CMI.
#  4. Compute CCA scores using CMI coefficients.
#  5. Evaluate brain–behavior correlation (Pearson + permutation).
#  6. Save predicted scores.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear all variables
rm(list = ls())

# Load required libraries
library(readxl)    # For reading Excel files
library(permute)   # For permutation utilities (shuffleSet)

# ------------------------------------------------------------------
# 2. Permutation Test Function
# ------------------------------------------------------------------

# Custom permutation test function to assess significance of correlation
my.perm <- function(numsub, nperm, X, Y) {
  actual_r <- as.numeric(cor(X, Y))             # Actual correlation
  perm.idx <- shuffleSet(numsub, nperm)         # Generate permutation indices
  r <- c()
  for (i in 1:nperm) {                          # Loop through permutations
    idx <- perm.idx[i, ]
    r <- c(r, cor(X[idx], Y))
  }
  p_perm <- length(which(r > actual_r)) / nperm # Proportion of permuted r > actual r
  p_perm
}
nperm <- 5000  # Number of permutations

# ------------------------------------------------------------------
# 3. Load Stanford Brain and Behavioral Data
# ------------------------------------------------------------------

# ROI labels from BrainNetome atlas
bn_atlas_file <- "data/atlas/bn_atlas.xlsx"
bn_atlas <- read_excel(bn_atlas_file, sheet = 1)
roi.names <- bn_atlas$Description[1:218]  # ROI names (first 218)

# Load Stanford behavioral data
beh_file <- "data/subjectlist/subjectlist_stanford_n231.csv"
beh <- read.csv(beh_file)

# Load Stanford brain GMV data
brain_file <- "data/gmv/gmv_stanford_n231.csv"
brain <- read.csv(brain_file)
brain <- as.matrix(brain[, -c(1:3)])  # Remove non-GMV columns (e.g., IDs)
sum(complete.cases(brain))            # Check for missing values
numsub <- dim(brain)[1]               # Number of subjects
numroi <- dim(brain)[2]               # Number of ROIs

# ------------------------------------------------------------------
# 4. Apply CMI Math CCA Model to Stanford Data
# ------------------------------------------------------------------

# Load original CMI math CCA model
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

# Project Stanford brain data into the PCA space from CMI
# (Subtract mean and multiply by PCA rotation matrix)
su_pcs <- sweep(brain, 2, my.pca$center, '-') %*% my.pca$rotation
su_pcs <- su_pcs[, 1:my.pca$numPC]  # Keep top PCs from CMI model

# Compute CCA mode 2 brain score
brain_score <- su_pcs %*% res$xcoef.raw[, 2]

# Compute behavioral score using CMI Y coefficients
beh_score <- beh$age * res$ycoef.raw[1, 2] +
  beh$mathrea_std * res$ycoef.raw[2, 2] +
  beh$numop_std   * res$ycoef.raw[3, 2]

# Correlation test (Pearson)
cor.test(brain_score, beh_score)

# Permutation test for significance
p.perm <- my.perm(numsub, nperm, brain_score, beh_score)
p.perm

# Save predicted brain and behavior scores
df <- data.frame(predicted_brain_score = brain_score, predicted_beh_score = beh_score)
write.csv(df, "results/cca/prediction/prediction_math.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 5. Apply CMI Reading CCA Model to Stanford Data
# ------------------------------------------------------------------

# Load original CMI reading CCA model
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

# Project Stanford brain data into the PCA space from CMI
su_pcs <- sweep(brain, 2, my.pca$center, '-') %*% my.pca$rotation
su_pcs <- su_pcs[, 1:my.pca$numPC]

# Compute CCA mode 2 brain score
brain_score <- su_pcs %*% res$xcoef.raw[, 2]

# Compute behavioral score using CMI Y coefficients
beh_score <- beh$age * res$ycoef.raw[1, 2] +
  beh$readcomp_std  * res$ycoef.raw[2, 2] +
  beh$wordread_std  * res$ycoef.raw[3, 2]

# Correlation test (Pearson)
cor.test(brain_score, beh_score)

# Permutation test for significance
p.perm <- my.perm(numsub, nperm, brain_score, beh_score)
p.perm

# Save predicted brain and behavior scores
df <- data.frame(predicted_brain_score = brain_score, predicted_beh_score = beh_score)
write.csv(df, "results/cca/prediction/prediction_reading.csv", row.names = FALSE)
