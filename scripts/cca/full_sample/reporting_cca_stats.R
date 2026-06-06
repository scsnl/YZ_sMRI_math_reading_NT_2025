############################################################
# CCA Reporting Full Statistics
# Author: Yuan Zhang
# Date: 2026-05-25
#
# Description:
# This script loads saved CCA RData files and generates
# reporting-ready statistics for editor requirements:
#   1. Wilks' Lambda dimension tests
#   2. Pillai's trace model fit tests
#   3. Canonical correlations with exact permutation p values
#      and bootstrap 95% confidence intervals
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

rm(list = ls())

library(CCA)
library(CCP)

out_dir <- "results/cca/reporting"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

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

format_num <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), x)
}

format_ci <- function(low, high, digits = 2) {
  sprintf(
    paste0("[%.", digits, "f, %.", digits, "f]"),
    low,
    high
  )
}

# ------------------------------------------------------------------
# 3. Bootstrap CI for canonical correlations
# ------------------------------------------------------------------

get_cancor_safe <- function(X, Y) {
  out <- tryCatch(
    {
      cc(X, Y)$cor
    },
    error = function(e) {
      rep(NA_real_, min(ncol(X), ncol(Y)))
    }
  )
  
  out
}

bootstrap_cancor_ci <- function(X, Y, nboot = 5000, conf = 0.95, seed = 1234) {
  set.seed(seed)
  
  n <- nrow(X)
  nmodes <- min(ncol(X), ncol(Y))
  
  boot_r <- matrix(NA_real_, nrow = nboot, ncol = nmodes)
  
  for (b in seq_len(nboot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    Xb <- X[idx, , drop = FALSE]
    Yb <- Y[idx, , drop = FALSE]
    
    r_b <- get_cancor_safe(Xb, Yb)
    
    if (length(r_b) >= nmodes) {
      boot_r[b, ] <- r_b[seq_len(nmodes)]
    }
  }
  
  alpha <- 1 - conf
  
  ci_low <- apply(
    boot_r,
    2,
    quantile,
    probs = alpha / 2,
    na.rm = TRUE
  )
  
  ci_high <- apply(
    boot_r,
    2,
    quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  
  data.frame(
    Mode = seq_len(nmodes),
    CI_low = ci_low,
    CI_high = ci_high,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 4. Extract Wilks' Lambda dimension tests
# ------------------------------------------------------------------

extract_wilks_tests <- function(cancor, n, p, q) {
  wilks <- p.asym(
    cancor,
    n,
    p,
    q,
    tstat = "Wilks"
  )
  
  wilks_df <- as.data.frame(wilks)
  colnames(wilks_df) <- make.names(colnames(wilks_df))
  
  # Print column names once if debugging is needed.
  # print(colnames(wilks_df))
  
  wilks_df$Dimension <- seq_len(nrow(wilks_df))
  
  stat_col <- grep("stat|Wilks|lambda", names(wilks_df), ignore.case = TRUE, value = TRUE)[1]
  f_col    <- grep("approx|F", names(wilks_df), ignore.case = TRUE, value = TRUE)[1]
  df1_col  <- grep("df1", names(wilks_df), ignore.case = TRUE, value = TRUE)[1]
  df2_col  <- grep("df2", names(wilks_df), ignore.case = TRUE, value = TRUE)[1]
  
  p_cols <- grep("p", names(wilks_df), ignore.case = TRUE, value = TRUE)
  p_col <- p_cols[length(p_cols)]
  
  if (any(is.na(c(stat_col, f_col, df1_col, df2_col, p_col)))) {
    stop(
      paste(
        "Could not identify all columns from CCP::p.asym output. Columns are:",
        paste(names(wilks_df), collapse = ", ")
      )
    )
  }
  
  out <- data.frame(
    Dimension = wilks_df$Dimension,
    Modes_tested = paste0(
      "Modes ",
      wilks_df$Dimension,
      " to ",
      nrow(wilks_df)
    ),
    Wilks_Lambda = wilks_df[[stat_col]],
    F = wilks_df[[f_col]],
    df1 = wilks_df[[df1_col]],
    df2 = wilks_df[[df2_col]],
    p = wilks_df[[p_col]],
    stringsAsFactors = FALSE
  )
  
  out$Full_report <- sprintf(
    "Wilks' Lambda = %.3f, F(%d, %d) = %.2f, p = %s",
    out$Wilks_Lambda,
    round(out$df1),
    round(out$df2),
    out$F,
    sapply(out$p, format_p)
  )
  
  out
}

# ------------------------------------------------------------------
# 5. Extract Pillai's trace model fit
# ------------------------------------------------------------------

extract_pillai_fit <- function(res, X, Y, check_saved = TRUE) {
  n <- nrow(X)
  pp <- ncol(Y)
  qq <- ncol(X)
  s <- min(pp, qq)
  df1_base <- max(pp, qq)
  df2_base <- n - max(pp, qq) - 1
  
  Fval <- (res$Pillai * df2_base) / ((s - res$Pillai) * df1_base)
  df1 <- s * df1_base
  df2 <- s * df2_base
  p_asym <- pf(Fval, df1, df2, lower.tail = FALSE)
  
  if (check_saved && !is.null(res$F.Pillai.perm)) {
    epsilon <- sqrt(.Machine$double.eps)
    
    p_recomputed <- (
      sum(res$F.Pillai.perm >= (Fval - epsilon), na.rm = TRUE) + 1
    ) / (length(res$F.Pillai.perm) + 1)
    
    if (abs(p_recomputed - res$pillai.p.perm) > 1e-10) {
      warning(
        "Saved Pillai permutation p value does not match recomputed +1 corrected p value."
      )
    }
  }
  
  data.frame(
    n = n,
    p_X = qq,
    p_Y = pp,
    Pillai_trace = res$Pillai,
    F = Fval,
    df1 = df1,
    df2 = df2,
    asymptotic_p = p_asym,
    permutation_p = res$pillai.p.perm,
    Full_report = sprintf(
      "Pillai's trace = %.3f, permutation p = %s",
      res$Pillai,
      format_p(res$pillai.p.perm)
    ),
    Full_report_with_asymptotic = sprintf(
      "Pillai's trace = %.3f, F(%d, %d) = %.2f, asymptotic p = %s, permutation p = %s",
      res$Pillai,
      round(df1),
      round(df2),
      Fval,
      format_p(p_asym),
      format_p(res$pillai.p.perm)
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 6. Extract canonical correlation permutation statistics
# ------------------------------------------------------------------

extract_cancor_perm_exact <- function(res, check_saved = TRUE) {
  observed <- res$cancor
  
  # Use saved permutation p values from updated mycca_noscale().
  p_exact <- res$cancor.p.perm
  p_fdr <- res$cancor.p.perm.adjust
  
  # Optional consistency check using saved null distribution.
  if (check_saved && !is.null(res$cancor.rand)) {
    null <- as.matrix(res$cancor.rand)
    nperm <- res$nperm
    
    p_recomputed <- sapply(seq_along(observed), function(i) {
      (sum(null[, i] >= observed[i], na.rm = TRUE) + 1) / (nperm + 1)
    })
    
    if (any(abs(p_recomputed - p_exact) > 1e-10, na.rm = TRUE)) {
      warning(
        "Saved canonical correlation permutation p values do not match recomputed +1 corrected p values."
      )
    }
  }
  
  data.frame(
    Mode = seq_along(observed),
    Canonical_r = observed,
    permutation_p = p_exact,
    FDR_corrected_permutation_p = p_fdr,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 7. Define saved CCA files
# ------------------------------------------------------------------

cca_files <- data.frame(
  Cohort = c("CMI-HBN", "CMI-HBN", "CMI-HBN", "Stanford", "Stanford", "Stanford"),
  Domain = c("Math", "Reading", "Joint", "Math", "Reading", "Joint"),
  File = c(
    "results/cca/cmi/wholebrain_cca_cmi_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    "results/cca/cmi/wholebrain_cca_cmi_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    "results/cca/stanford/wholebrain_cca_stanford_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    "results/cca/stanford/wholebrain_cca_stanford_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData",
    "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined/CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# 8. Loop over files and extract statistics
# ------------------------------------------------------------------

all_wilks <- list()
all_pillai <- list()
all_cancor <- list()

nboot_ci <- 5000

for (i in seq_len(nrow(cca_files))) {
  cohort <- cca_files$Cohort[i]
  domain <- cca_files$Domain[i]
  file <- cca_files$File[i]
  
  cat("\n----------------------------------------\n")
  cat("Processing:", cohort, domain, "\n")
  cat("File:", file, "\n")
  cat("----------------------------------------\n")
  
  if (!file.exists(file)) {
    stop(paste("File not found:", file))
  }
  
  # Loads saved objects: data, Xorig, X, Y, res, my.pca
  load(file)
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  # --------------------------
  # Wilks' Lambda dimension tests
  # --------------------------
  
  wilks_df <- extract_wilks_tests(
    cancor = res$cancor,
    n = n,
    p = p,
    q = q
  )
  
  wilks_df$Cohort <- cohort
  wilks_df$Domain <- domain
  wilks_df$n <- n
  wilks_df$p_X <- p
  wilks_df$p_Y <- q
  
  wilks_df <- wilks_df[, c(
    "Cohort",
    "Domain",
    "n",
    "p_X",
    "p_Y",
    "Dimension",
    "Modes_tested",
    "Wilks_Lambda",
    "F",
    "df1",
    "df2",
    "p",
    "Full_report"
  )]
  
  all_wilks[[paste(cohort, domain, sep = "_")]] <- wilks_df
  
  # --------------------------
  # Pillai's trace model fit
  # --------------------------
  
  pillai_df <- extract_pillai_fit(
    res = res,
    X = X,
    Y = Y,
    check_saved = TRUE
  )
  
  pillai_df$Cohort <- cohort
  pillai_df$Domain <- domain
  
  pillai_df <- pillai_df[, c(
    "Cohort",
    "Domain",
    "n",
    "p_X",
    "p_Y",
    "Pillai_trace",
    "F",
    "df1",
    "df2",
    "asymptotic_p",
    "permutation_p",
    "Full_report",
    "Full_report_with_asymptotic"
  )]
  
  all_pillai[[paste(cohort, domain, sep = "_")]] <- pillai_df
  
  # --------------------------
  # Canonical correlations
  # --------------------------
  
  cancor_df <- extract_cancor_perm_exact(
    res = res,
    check_saved = TRUE
  )
  
  cat("Bootstrapping canonical correlation 95% CIs...\n")
  
  ci_df <- bootstrap_cancor_ci(
    X = X,
    Y = Y,
    nboot = nboot_ci,
    conf = 0.95,
    seed = 1000 + i
  )
  
  cancor_df <- merge(
    cancor_df,
    ci_df,
    by = "Mode",
    all.x = TRUE
  )
  
  cancor_df$Cohort <- cohort
  cancor_df$Domain <- domain
  cancor_df$n <- n
  cancor_df$p_X <- p
  cancor_df$p_Y <- q
  
  cancor_df$Full_report <- sprintf(
    "canonical r = %.2f, 95%% CI %s, permutation p = %s, FDR-corrected permutation p = %s",
    cancor_df$Canonical_r,
    format_ci(cancor_df$CI_low, cancor_df$CI_high),
    sapply(cancor_df$permutation_p, format_p),
    sapply(cancor_df$FDR_corrected_permutation_p, format_p)
  )
  
  cancor_df <- cancor_df[, c(
    "Cohort",
    "Domain",
    "n",
    "p_X",
    "p_Y",
    "Mode",
    "Canonical_r",
    "CI_low",
    "CI_high",
    "permutation_p",
    "FDR_corrected_permutation_p",
    "Full_report"
  )]
  
  all_cancor[[paste(cohort, domain, sep = "_")]] <- cancor_df
}

# ------------------------------------------------------------------
# 9. Combine and save outputs
# ------------------------------------------------------------------

cca_wilks_dimension_tests <- do.call(rbind, all_wilks)
cca_pillai_model_fit <- do.call(rbind, all_pillai)
cca_canonical_correlations <- do.call(rbind, all_cancor)

write.csv(
  cca_wilks_dimension_tests,
  file.path(out_dir, "cca_wilks_dimension_tests_full_statistics.csv"),
  row.names = FALSE
)

write.csv(
  cca_pillai_model_fit,
  file.path(out_dir, "cca_pillai_model_fit_full_statistics.csv"),
  row.names = FALSE
)

write.csv(
  cca_canonical_correlations,
  file.path(out_dir, "cca_canonical_correlation_full_statistics.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# 10. Print results
# ------------------------------------------------------------------

cat("\n### Wilks' Lambda Dimension Tests ###\n")
print(cca_wilks_dimension_tests)

cat("\n### Pillai's Trace Model Fit Tests ###\n")
print(cca_pillai_model_fit)

cat("\n### Canonical Correlations ###\n")
print(cca_canonical_correlations)

cat("\nCCA reporting statistics completed.\n")