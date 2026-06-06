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

# ------------------------------------------------------------------
# 2. Load original CCA results and subjectlists
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
# 3. Load subjectlists with SST_GVOL and SST_BVOL
# ------------------------------------------------------------------
cmi_subj  <- read.csv("data/subjectlist/subjectlist_cmi_n760.csv")
stan_subj <- read.csv("data/subjectlist/subjectlist_stanford_n231.csv")

# Make Stanford IDs numeric if needed
stan_subj$pid <- as.numeric(stan_subj$pid)
stan_math_id  <- as.numeric(stan_math_id)
stan_read_id  <- as.numeric(stan_read_id)

# ------------------------------------------------------------------
# 4. Helper function: permutation test for partial correlation
# ------------------------------------------------------------------
partial_cor_perm_test <- function(U, V, COVAR, nperm = 5000) {
  obs <- ppcor::pcor.test(U, V, COVAR, method = "pearson")
  obs_r <- as.numeric(obs$estimate)
  
  perm_r <- numeric(nperm)
  for (i in seq_len(nperm)) {
    V_perm <- sample(V, length(V), replace = FALSE)
    perm_r[i] <- as.numeric(ppcor::pcor.test(U, V_perm, COVAR, method = "pearson")$estimate)
  }
  
  p_perm <- (sum(abs(perm_r) >= abs(obs_r)) + 1) / (nperm + 1)
  
  list(
    obs_r = obs_r,
    asym_p = obs$p.value,
    perm_p = p_perm,
    perm_r = perm_r
  )
}

# ------------------------------------------------------------------
# 5. Main wrapper
# ------------------------------------------------------------------
run_partial_corr <- function(ids, U, V, subj_df, id_col, covar_col,
                             dataset_name, domain_name, nperm = 5000) {
  
  # Match IDs
  idx <- match(ids, subj_df[[id_col]])
  covar <- subj_df[[covar_col]][idx]
  
  # Keep complete cases
  keep <- complete.cases(U, V, covar)
  U2 <- U[keep]
  V2 <- V[keep]
  covar2 <- covar[keep]
  
  # Pearson correlation
  pearson_res <- cor.test(U2, V2, method = "pearson")
  
  # Partial correlation
  pcor_res <- ppcor::pcor.test(U2, V2, covar2, method = "pearson")
  
  # Permutation test for partial correlation
  perm_res <- partial_cor_perm_test(U2, V2, covar2, nperm = nperm)
  
  # Save permutation distribution
  perm_df <- data.frame(perm_partial_r = perm_res$perm_r)
  out_perm_file <- paste0(
    "results/cca/partial_corr/permutation_",
    gsub("-", "", dataset_name), "_",
    tolower(domain_name), "_partialcorr_", covar_col, ".csv"
  )
  write.csv(perm_df, out_perm_file, row.names = FALSE)
  
  data.frame(
    Dataset = dataset_name,
    Domain = domain_name,
    Covariate = covar_col,
    N = length(U2),
    Pearson_r = unname(pearson_res$estimate),
    Pearson_p = pearson_res$p.value,
    Partial_r = as.numeric(pcor_res$estimate),
    Partial_p = pcor_res$p.value,
    Partial_perm_p = perm_res$perm_p,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 6. Run analyses: SST_GVOL
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
  nperm = 5000
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
  nperm = 5000
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
  nperm = 5000
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
  nperm = 5000
)

# ------------------------------------------------------------------
# 7. Run analyses: SST_BVOL
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
  nperm = 5000
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
  nperm = 5000
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
  nperm = 5000
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
  nperm = 5000
)

all_results <- rbind(
  res_cmi_math_gvol, res_cmi_read_gvol, res_stan_math_gvol, res_stan_read_gvol,
  res_cmi_math_bvol, res_cmi_read_bvol, res_stan_math_bvol, res_stan_read_bvol
)

# ------------------------------------------------------------------
# 8. Print and save results
# ------------------------------------------------------------------
print(all_results)

write.csv(
  all_results,
  "results/cca/partial_corr/original_mode_partialcorr_control_SST_GVOL_SST_BVOL_with_perm.csv",
  row.names = FALSE
)

cat("\nSaved results to: results/cca/original_mode_partialcorr_control_SST_GVOL_SST_BVOL_with_perm.csv\n")