############################################################
# Partial correlation between original U and V controlling
# for global morphometric covariates
#   1) SST_GVOL
#   2) SST_BVOL
# + permutation test for significance
#
# Author: Yuan Zhang
# Date: 2026-03-30
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

library(ppcor)

set.seed(1234)

dir.create("results/cca/partial_corr", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------

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

# Bootstrap CI for partial correlation
bootstrap_partial_ci <- function(U, V, COVAR, nboot = 5000, conf = 0.95, seed = 1234) {
  set.seed(seed)
  
  keep <- complete.cases(U, V, COVAR)
  U <- as.numeric(U[keep])
  V <- as.numeric(V[keep])
  COVAR <- as.numeric(COVAR[keep])
  
  n <- length(U)
  boot_r <- numeric(nboot)
  
  for (b in seq_len(nboot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    boot_r[b] <- as.numeric(
      ppcor::pcor.test(U[idx], V[idx], COVAR[idx], method = "pearson")$estimate
    )
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

# ------------------------------------------------------------------
# 3. Load original CCA results and subjectlists
# ------------------------------------------------------------------

###############################
## CMI - Math original model ##
###############################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

cmi_math_id <- data$beh$oak_id
cmi_math_U  <- res$Cx[, 2]
cmi_math_V  <- res$Cy[, 2]

#################################
## CMI - Reading original model ##
#################################
orig_path <- "results/cca/cmi/wholebrain_cca_cmi_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

cmi_read_id <- data$beh$oak_id
cmi_read_U  <- res$Cx[, 2]
cmi_read_V  <- res$Cy[, 2]

####################################
## Stanford - Math original model ##
####################################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_math/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

stan_math_id <- data$beh$pid
stan_math_U  <- res$Cx[, 2]
stan_math_V  <- res$Cy[, 2]

######################################
## Stanford - Reading original model ##
######################################
orig_path <- "results/cca/stanford/wholebrain_cca_stanford_reading/"
orig_file <- "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste0(orig_path, orig_file))

stan_read_id <- data$beh$pid
stan_read_U  <- res$Cx[, 2]
stan_read_V  <- res$Cy[, 2]

# ------------------------------------------------------------------
# 4. Load subjectlists with SST_GVOL and SST_BVOL
# ------------------------------------------------------------------

cmi_subj  <- read.csv("data/subjectlist/subjectlist_cmi_n760.csv")
stan_subj <- read.csv("data/subjectlist/subjectlist_stanford_n231.csv")

stan_subj$pid <- as.numeric(stan_subj$pid)
stan_math_id  <- as.numeric(stan_math_id)
stan_read_id  <- as.numeric(stan_read_id)

# ------------------------------------------------------------------
# 5. Helper function: permutation test for partial correlation
# ------------------------------------------------------------------

partial_cor_perm_test <- function(U, V, COVAR, nperm = 5000) {
  keep <- complete.cases(U, V, COVAR)
  U <- as.numeric(U[keep])
  V <- as.numeric(V[keep])
  COVAR <- as.numeric(COVAR[keep])
  
  obs <- ppcor::pcor.test(U, V, COVAR, method = "pearson")
  obs_r <- as.numeric(obs$estimate)
  
  perm_r <- numeric(nperm)
  
  for (i in seq_len(nperm)) {
    V_perm <- sample(V, length(V), replace = FALSE)
    perm_r[i] <- as.numeric(
      ppcor::pcor.test(U, V_perm, COVAR, method = "pearson")$estimate
    )
  }
  
  # Two-sided permutation p with +1 correction
  p_perm <- (sum(abs(perm_r) >= abs(obs_r)) + 1) / (nperm + 1)
  
  list(
    obs_r = obs_r,
    asym_p = obs$p.value,
    perm_p = p_perm,
    perm_r = perm_r
  )
}

# ------------------------------------------------------------------
# 6. Main wrapper
# ------------------------------------------------------------------

run_partial_corr <- function(ids, U, V, subj_df, id_col, covar_col,
                             dataset_name, domain_name, nperm = 5000,
                             nboot = 5000) {
  
  # Match IDs
  idx <- match(ids, subj_df[[id_col]])
  covar <- subj_df[[covar_col]][idx]
  
  # Keep complete cases
  keep <- complete.cases(U, V, covar)
  U2 <- as.numeric(U[keep])
  V2 <- as.numeric(V[keep])
  covar2 <- as.numeric(covar[keep])
  
  n <- length(U2)
  
  # ------------------------------------------------------------
  # Pearson correlation between original U and V
  # ------------------------------------------------------------
  pearson_res <- cor.test(U2, V2, method = "pearson")
  
  pearson_t <- unname(pearson_res$statistic)
  pearson_df <- unname(pearson_res$parameter)
  pearson_r <- unname(pearson_res$estimate)
  pearson_p <- pearson_res$p.value
  pearson_ci_low <- pearson_res$conf.int[1]
  pearson_ci_high <- pearson_res$conf.int[2]
  
  # ------------------------------------------------------------
  # Partial correlation controlling for covariate
  # ------------------------------------------------------------
  pcor_res <- ppcor::pcor.test(U2, V2, covar2, method = "pearson")
  
  partial_r <- as.numeric(pcor_res$estimate)
  partial_p <- pcor_res$p.value
  
  # For one covariate, df = n - number_of_covariates - 2 = n - 3
  n_covar <- 1
  partial_df <- n - n_covar - 2
  partial_t <- partial_r * sqrt(partial_df / (1 - partial_r^2))
  
  # Bootstrap 95% CI for partial r
  boot_ci <- bootstrap_partial_ci(
    U = U2,
    V = V2,
    COVAR = covar2,
    nboot = nboot,
    conf = 0.95,
    seed = 1000 + nchar(dataset_name) + nchar(domain_name) + nchar(covar_col)
  )
  
  # ------------------------------------------------------------
  # Permutation test for partial correlation
  # ------------------------------------------------------------
  perm_res <- partial_cor_perm_test(U2, V2, covar2, nperm = nperm)
  
  # Save permutation distribution
  perm_df <- data.frame(perm_partial_r = perm_res$perm_r)
  
  out_perm_file <- paste0(
    "results/cca/partial_corr/permutation_",
    gsub("-", "", dataset_name), "_",
    tolower(domain_name), "_partialcorr_", covar_col, ".csv"
  )
  
  write.csv(perm_df, out_perm_file, row.names = FALSE)
  
  # Save bootstrap distribution
  boot_df <- data.frame(boot_partial_r = boot_ci$boot_r)
  
  out_boot_file <- paste0(
    "results/cca/partial_corr/bootstrap_",
    gsub("-", "", dataset_name), "_",
    tolower(domain_name), "_partialcorr_", covar_col, ".csv"
  )
  
  write.csv(boot_df, out_boot_file, row.names = FALSE)
  
  # ------------------------------------------------------------
  # Return reporting-ready result
  # ------------------------------------------------------------
  data.frame(
    Dataset = dataset_name,
    Domain = domain_name,
    Covariate = covar_col,
    N = n,
    
    Pearson_t = pearson_t,
    Pearson_df = pearson_df,
    Pearson_r = pearson_r,
    Pearson_CI_low = pearson_ci_low,
    Pearson_CI_high = pearson_ci_high,
    Pearson_p = pearson_p,
    
    Partial_t = partial_t,
    Partial_df = partial_df,
    Partial_r = partial_r,
    Partial_CI_low = boot_ci$ci_low,
    Partial_CI_high = boot_ci$ci_high,
    Partial_asymptotic_p = partial_p,
    Partial_perm_p = perm_res$perm_p,
    
    Full_report_pearson = sprintf(
      "t(%d) = %.2f, p = %s, r = %.2f, 95%% CI %s",
      round(pearson_df),
      pearson_t,
      format_p(pearson_p),
      pearson_r,
      format_ci(pearson_ci_low, pearson_ci_high)
    ),
    
    Full_report_partial_asymptotic = sprintf(
      "t(%d) = %.2f, p = %s, partial r = %.2f, bootstrap 95%% CI %s",
      round(partial_df),
      partial_t,
      format_p(partial_p),
      partial_r,
      format_ci(boot_ci$ci_low, boot_ci$ci_high)
    ),
    
    Full_report_partial_perm = sprintf(
      "t(%d) = %.2f, partial r = %.2f, bootstrap 95%% CI %s, two-sided permutation p = %s",
      round(partial_df),
      partial_t,
      partial_r,
      format_ci(boot_ci$ci_low, boot_ci$ci_high),
      format_p(perm_res$perm_p)
    ),
    
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 7. Run analyses: SST_GVOL
# ------------------------------------------------------------------

res_cmi_math_gvol <- run_partial_corr(
  ids = cmi_math_id,
  U = cmi_math_U,
  V = cmi_math_V,
  subj_df = cmi_subj,
  id_col = "oak_id",
  covar_col = "SST_GVOL",
  dataset_name = "CMI-HBN",
  domain_name = "Math",
  nperm = 5000,
  nboot = 5000
)

res_cmi_read_gvol <- run_partial_corr(
  ids = cmi_read_id,
  U = cmi_read_U,
  V = cmi_read_V,
  subj_df = cmi_subj,
  id_col = "oak_id",
  covar_col = "SST_GVOL",
  dataset_name = "CMI-HBN",
  domain_name = "Reading",
  nperm = 5000,
  nboot = 5000
)

res_stan_math_gvol <- run_partial_corr(
  ids = stan_math_id,
  U = stan_math_U,
  V = stan_math_V,
  subj_df = stan_subj,
  id_col = "pid",
  covar_col = "SST_GVOL",
  dataset_name = "Stanford",
  domain_name = "Math",
  nperm = 5000,
  nboot = 5000
)

res_stan_read_gvol <- run_partial_corr(
  ids = stan_read_id,
  U = stan_read_U,
  V = stan_read_V,
  subj_df = stan_subj,
  id_col = "pid",
  covar_col = "SST_GVOL",
  dataset_name = "Stanford",
  domain_name = "Reading",
  nperm = 5000,
  nboot = 5000
)

# ------------------------------------------------------------------
# 8. Run analyses: SST_BVOL
# ------------------------------------------------------------------

res_cmi_math_bvol <- run_partial_corr(
  ids = cmi_math_id,
  U = cmi_math_U,
  V = cmi_math_V,
  subj_df = cmi_subj,
  id_col = "oak_id",
  covar_col = "SST_BVOL",
  dataset_name = "CMI-HBN",
  domain_name = "Math",
  nperm = 5000,
  nboot = 5000
)

res_cmi_read_bvol <- run_partial_corr(
  ids = cmi_read_id,
  U = cmi_read_U,
  V = cmi_read_V,
  subj_df = cmi_subj,
  id_col = "oak_id",
  covar_col = "SST_BVOL",
  dataset_name = "CMI-HBN",
  domain_name = "Reading",
  nperm = 5000,
  nboot = 5000
)

res_stan_math_bvol <- run_partial_corr(
  ids = stan_math_id,
  U = stan_math_U,
  V = stan_math_V,
  subj_df = stan_subj,
  id_col = "pid",
  covar_col = "SST_BVOL",
  dataset_name = "Stanford",
  domain_name = "Math",
  nperm = 5000,
  nboot = 5000
)

res_stan_read_bvol <- run_partial_corr(
  ids = stan_read_id,
  U = stan_read_U,
  V = stan_read_V,
  subj_df = stan_subj,
  id_col = "pid",
  covar_col = "SST_BVOL",
  dataset_name = "Stanford",
  domain_name = "Reading",
  nperm = 5000,
  nboot = 5000
)

all_results <- rbind(
  res_cmi_math_gvol,
  res_cmi_read_gvol,
  res_stan_math_gvol,
  res_stan_read_gvol,
  res_cmi_math_bvol,
  res_cmi_read_bvol,
  res_stan_math_bvol,
  res_stan_read_bvol
)

# ------------------------------------------------------------------
# 9. Print and save results
# ------------------------------------------------------------------

print(all_results)

write.csv(
  all_results,
  "results/cca/partial_corr/original_mode_partialcorr_control_SST_GVOL_SST_BVOL_with_perm_full_statistics.csv",
  row.names = FALSE
)

cat("\nSaved results to: results/cca/partial_corr/original_mode_partialcorr_control_SST_GVOL_SST_BVOL_with_perm_full_statistics.csv\n")