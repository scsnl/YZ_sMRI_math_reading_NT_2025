############################################################
# Prediction Using CCA Models (CMI → Stanford)
# Author: Yuan Zhang
# Date: 2026-05-25
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
#  5. Evaluate brain–behavior correlation using Pearson r,
#     bootstrap 95% CI, and permutation testing.
#  6. Save predicted scores and reporting-ready statistics.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

rm(list = ls())

library(readxl)
library(permute)

# ------------------------------------------------------------------
# 2. Custom functions
# ------------------------------------------------------------------

my.perm <- function(numsub, nperm, X, Y) {
  keep <- complete.cases(X, Y)
  X <- as.numeric(X[keep])
  Y <- as.numeric(Y[keep])
  numsub <- length(X)
  
  actual_r <- as.numeric(cor(X, Y))
  perm.idx <- shuffleSet(numsub, nperm)
  
  r <- c()
  for (i in 1:nperm) {
    idx <- perm.idx[i, ]
    r <- c(r, cor(X[idx], Y))
  }
  
  # Upper-tail permutation p with +1 correction
  p_perm <- (sum(r >= actual_r) + 1) / (nperm + 1)
  
  list(
    actual_r = actual_r,
    p_perm = p_perm,
    perm_r = r
  )
}

bootstrap_cor_ci <- function(X, Y, nboot = 5000, conf = 0.95, seed = 1234) {
  set.seed(seed)
  
  keep <- complete.cases(X, Y)
  X <- as.numeric(X[keep])
  Y <- as.numeric(Y[keep])
  
  n <- length(X)
  boot_r <- numeric(nboot)
  
  for (b in seq_len(nboot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    boot_r[b] <- cor(X[idx], Y[idx])
  }
  
  alpha <- 1 - conf
  
  ci <- quantile(
    boot_r,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE
  )
  
  list(
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    boot_r = boot_r
  )
}

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

extract_prediction_stats <- function(ct, p.perm, boot.ci, domain, mode = 2) {
  data.frame(
    Domain = domain,
    Mode = mode,
    n = unname(ct$parameter) + 2,
    Statistic = sprintf(
      "t(%d) = %.2f",
      round(unname(ct$parameter)),
      unname(ct$statistic)
    ),
    Effect_size = sprintf("r = %.2f", unname(ct$estimate)),
    CI_95 = sprintf(
      "bootstrap 95%% CI %s",
      format_ci(boot.ci$ci_low, boot.ci$ci_high)
    ),
    Parametric_p = ct$p.value,
    Permutation_p = p.perm$p_perm,
    Full_report = sprintf(
      "t(%d) = %.2f, r = %.2f, bootstrap 95%% CI %s, permutation p = %s",
      round(unname(ct$parameter)),
      unname(ct$statistic),
      unname(ct$estimate),
      format_ci(boot.ci$ci_low, boot.ci$ci_high),
      format_p(p.perm$p_perm)
    )
  )
}

# ------------------------------------------------------------------
# 3. Load Stanford Brain and Behavioral Data
# ------------------------------------------------------------------

bn_atlas_file <- "data/atlas/bn_atlas.xlsx"
bn_atlas <- read_excel(bn_atlas_file, sheet = 1)
roi.names <- bn_atlas$Description[1:218]

beh_file <- "data/subjectlist/subjectlist_stanford_n231.csv"
beh <- read.csv(beh_file)

brain_file <- "data/gmv/gmv_stanford_n231.csv"
brain <- read.csv(brain_file)
brain <- as.matrix(brain[, -c(1:3)])
colnames(brain) <- roi.names

sum(complete.cases(brain))
numsub <- dim(brain)[1]
numroi <- dim(brain)[2]

nperm <- 5000
nboot <- 5000

prediction_stats_list <- list()

# ------------------------------------------------------------------
# 4. Apply CMI Math CCA Model to Stanford Data
# ------------------------------------------------------------------

orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

su_pcs <- sweep(brain, 2, my.pca$center, "-") %*% my.pca$rotation
su_pcs <- su_pcs[, 1:my.pca$numPC]

brain_score <- su_pcs %*% res$xcoef.raw[, 2]

beh_score <- beh$age * res$ycoef.raw[1, 2] +
  beh$mathrea_std * res$ycoef.raw[2, 2] +
  beh$numop_std   * res$ycoef.raw[3, 2]

# Pearson correlation
ct_math <- cor.test(brain_score, beh_score)

# Upper-tail permutation test
p.perm_math <- my.perm(numsub, nperm, brain_score, beh_score)

# Bootstrap CI for Pearson r
boot.ci_math <- bootstrap_cor_ci(
  X = brain_score,
  Y = beh_score,
  nboot = nboot,
  conf = 0.95,
  seed = 1001
)

prediction_stats_list[["math"]] <- extract_prediction_stats(
  ct = ct_math,
  p.perm = p.perm_math,
  boot.ci = boot.ci_math,
  domain = "Math",
  mode = 2
)

print(ct_math)
print(p.perm_math$p_perm)
print(prediction_stats_list[["math"]]$Full_report)

df <- data.frame(
  predicted_brain_score = -1 * as.numeric(brain_score),
  predicted_beh_score = -1 * as.numeric(beh_score)
)

write.csv(
  df,
  "results/cca/prediction/prediction_math.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------
# 5. Apply CMI Reading CCA Model to Stanford Data
# ------------------------------------------------------------------

orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

su_pcs <- sweep(brain, 2, my.pca$center, "-") %*% my.pca$rotation
su_pcs <- su_pcs[, 1:my.pca$numPC]

brain_score <- su_pcs %*% res$xcoef.raw[, 2]

beh_score <- beh$age * res$ycoef.raw[1, 2] +
  beh$readcomp_std  * res$ycoef.raw[2, 2] +
  beh$wordread_std  * res$ycoef.raw[3, 2]

# Pearson correlation
ct_reading <- cor.test(brain_score, beh_score)

# Upper-tail permutation test
p.perm_reading <- my.perm(numsub, nperm, brain_score, beh_score)

# Bootstrap CI for Pearson r
boot.ci_reading <- bootstrap_cor_ci(
  X = brain_score,
  Y = beh_score,
  nboot = nboot,
  conf = 0.95,
  seed = 1002
)

prediction_stats_list[["reading"]] <- extract_prediction_stats(
  ct = ct_reading,
  p.perm = p.perm_reading,
  boot.ci = boot.ci_reading,
  domain = "Reading",
  mode = 2
)

print(ct_reading)
print(p.perm_reading$p_perm)
print(prediction_stats_list[["reading"]]$Full_report)

df <- data.frame(
  predicted_brain_score = -1 * as.numeric(brain_score),
  predicted_beh_score = -1 * as.numeric(beh_score)
)

write.csv(
  df,
  "results/cca/prediction/prediction_reading.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------
# 6. Summary table
# ------------------------------------------------------------------

prediction_stats <- do.call(rbind, prediction_stats_list)

# FDR correction across math and reading prediction tests
prediction_stats$FDR_corrected_permutation_p <- p.adjust(
  prediction_stats$Permutation_p,
  method = "fdr"
)

prediction_stats$Full_report_FDR <- sprintf(
  "%s, FDR-corrected permutation p = %s",
  prediction_stats$Full_report,
  sapply(prediction_stats$FDR_corrected_permutation_p, format_p)
)

write.csv(
  prediction_stats,
  "results/cca/prediction/prediction_statistics_full_report.csv",
  row.names = FALSE
)

print(prediction_stats)