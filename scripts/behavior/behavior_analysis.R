############################################################
# Behavioral Analysis Script
# Author: Yuan Zhang
# Date: 2026-05-25
#
# Description:
# This script performs behavioral analyses for two cohorts
# (CMI-HBN and Stanford), including descriptive statistics,
# Pearson correlations, and group comparisons. Inferential
# statistics are formatted to include statistic(df), p value,
# effect size, and 95% confidence intervals.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

rm(list = ls())

library(effectsize)

# Create output folder if needed
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
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

format_ci <- function(low, high, digits = 2) {
  sprintf(
    paste0("[%.", digits, "f, %.", digits, "f]"),
    low, high
  )
}

mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  sprintf(
    paste0("%.", digits, "f \u00B1 %.", digits, "f"),
    mean(x),
    sd(x)
  )
}

range_string <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  sprintf(
    paste0("%.", digits, "f\u2013%.", digits, "f"),
    min(x),
    max(x)
  )
}

get_descriptive_stats <- function(data, cohort_name, variables) {
  out <- do.call(
    rbind,
    lapply(names(variables), function(label) {
      var <- variables[[label]]
      x <- data[[var]]
      
      data.frame(
        Cohort = cohort_name,
        Variable = label,
        n = sum(!is.na(x)),
        Mean_SD = mean_sd(x),
        Range = range_string(x),
        stringsAsFactors = FALSE
      )
    })
  )
  
  out
}

get_pearson_stats <- function(data, cohort_name, xvar, yvar, x_label, y_label) {
  df_data <- data[, c(xvar, yvar)]
  df_data <- df_data[complete.cases(df_data), ]
  
  ct <- cor.test(
    df_data[[xvar]],
    df_data[[yvar]],
    method = "pearson"
  )
  
  stat_df <- round(unname(ct$parameter))
  stat_value <- unname(ct$statistic)
  r_value <- unname(ct$estimate)
  p_value <- ct$p.value
  ci_low <- ct$conf.int[1]
  ci_high <- ct$conf.int[2]
  
  data.frame(
    Cohort = cohort_name,
    Variable_1 = x_label,
    Variable_2 = y_label,
    n = nrow(df_data),
    Statistic = sprintf("t(%d) = %.2f", stat_df, stat_value),
    p = sprintf("p = %s", format_p(p_value)),
    Effect_size = sprintf("r = %.2f", r_value),
    CI_95 = sprintf("95%% CI %s", format_ci(ci_low, ci_high)),
    Full_report = sprintf(
      "t(%d) = %.2f, p = %s, r = %.2f, 95%% CI %s",
      stat_df,
      stat_value,
      format_p(p_value),
      r_value,
      format_ci(ci_low, ci_high)
    ),
    stringsAsFactors = FALSE
  )
}

get_cont_group_stats <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  # Equal variance t-test
  tt <- t.test(x, y, var.equal = TRUE)
  
  # Cohen's d using pooled SD, matching equal variance t-test.
  d <- cohens_d(x, y, pooled_sd = TRUE, ci = 0.95)
  
  data.frame(
    CMI = mean_sd(x),
    Stanford = mean_sd(y),
    Group_difference = sprintf(
      "t(%d) = %.2f, p = %s, Cohen's d = %.2f, 95%% CI %s",
      round(unname(tt$parameter)),
      unname(tt$statistic),
      format_p(tt$p.value),
      d$Cohens_d,
      format_ci(d$CI_low, d$CI_high)
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------

beh_cmi <- read.csv("data/subjectlist/subjectlist_cmi_n760.csv")
beh_su <- read.csv("data/subjectlist/subjectlist_stanford_n231.csv")

# Harmonize sex coding for group comparison
beh_cmi$genderF <- ifelse(
  beh_cmi$Sex == "F", 1,
  ifelse(beh_cmi$Sex == "M", 0, NA)
)

# ------------------------------------------------------------------
# 4. Descriptive statistics
# ------------------------------------------------------------------

cmi_desc_vars <- list(
  Age = "age",
  FSIQ = "WISC_FSIQ",
  Numerical_operations = "numop_std",
  Math_problem_solving = "mathprob_std",
  Word_reading = "wordread_std",
  Reading_comprehension = "readcomp_std"
)

stanford_desc_vars <- list(
  Age = "age",
  FSIQ = "fsiq",
  Numerical_operations = "numop_std",
  Math_reasoning = "mathrea_std",
  Word_reading = "wordread_std",
  Reading_comprehension = "readcomp_std"
)

cmi_desc <- get_descriptive_stats(
  data = beh_cmi,
  cohort_name = "CMI-HBN",
  variables = cmi_desc_vars
)

stanford_desc <- get_descriptive_stats(
  data = beh_su,
  cohort_name = "Stanford",
  variables = stanford_desc_vars
)

behavioral_descriptives <- rbind(cmi_desc, stanford_desc)

# Sex counts
sex_counts <- data.frame(
  Cohort = c("CMI-HBN", "Stanford"),
  Male = c(
    sum(beh_cmi$genderF == 0, na.rm = TRUE),
    sum(beh_su$genderF == 0, na.rm = TRUE)
  ),
  Female = c(
    sum(beh_cmi$genderF == 1, na.rm = TRUE),
    sum(beh_su$genderF == 1, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  behavioral_descriptives,
  "results/behavior/behavioral_descriptive_statistics.csv",
  row.names = FALSE
)

write.csv(
  sex_counts,
  "results/behavior/behavioral_sex_counts.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------
# 5. Pearson correlation analyses
# ------------------------------------------------------------------

cmi_cor_tests <- list(
  # Math-reading correlations
  list("numop_std", "wordread_std", "Numerical operations", "Word reading"),
  list("numop_std", "readcomp_std", "Numerical operations", "Reading comprehension"),
  list("mathprob_std", "wordread_std", "Math problem solving", "Word reading"),
  list("mathprob_std", "readcomp_std", "Math problem solving", "Reading comprehension"),
  
  # Age correlations
  list("age", "numop_std", "Age", "Numerical operations"),
  list("age", "mathprob_std", "Age", "Math problem solving"),
  list("age", "wordread_std", "Age", "Word reading"),
  list("age", "readcomp_std", "Age", "Reading comprehension"),
  
  # FSIQ correlations
  list("numop_std", "WISC_FSIQ", "Numerical operations", "FSIQ"),
  list("mathprob_std", "WISC_FSIQ", "Math problem solving", "FSIQ"),
  list("wordread_std", "WISC_FSIQ", "Word reading", "FSIQ"),
  list("readcomp_std", "WISC_FSIQ", "Reading comprehension", "FSIQ"),
  
  # Working memory correlations
  list("numop_std", "WISC_WMI", "Numerical operations", "Working memory index"),
  list("mathprob_std", "WISC_WMI", "Math problem solving", "Working memory index"),
  list("wordread_std", "WISC_WMI", "Word reading", "Working memory index"),
  list("readcomp_std", "WISC_WMI", "Reading comprehension", "Working memory index")
)

stanford_cor_tests <- list(
  # Math-reading correlations
  list("numop_std", "wordread_std", "Numerical operations", "Word reading"),
  list("numop_std", "readcomp_std", "Numerical operations", "Reading comprehension"),
  list("mathrea_std", "wordread_std", "Math reasoning", "Word reading"),
  list("mathrea_std", "readcomp_std", "Math reasoning", "Reading comprehension"),
  
  # Age correlations
  list("age", "numop_std", "Age", "Numerical operations"),
  list("age", "mathrea_std", "Age", "Math reasoning"),
  list("age", "wordread_std", "Age", "Word reading"),
  list("age", "readcomp_std", "Age", "Reading comprehension"),
  
  # FSIQ correlations
  list("numop_std", "fsiq", "Numerical operations", "FSIQ"),
  list("mathrea_std", "fsiq", "Math reasoning", "FSIQ"),
  list("wordread_std", "fsiq", "Word reading", "FSIQ"),
  list("readcomp_std", "fsiq", "Reading comprehension", "FSIQ"),
  
  # Working memory correlations
  list("numop_std", "Backward_digit_recall_std", "Numerical operations", "Backward digit recall"),
  list("mathrea_std", "Backward_digit_recall_std", "Math reasoning", "Backward digit recall"),
  list("wordread_std", "Backward_digit_recall_std", "Word reading", "Backward digit recall"),
  list("readcomp_std", "Backward_digit_recall_std", "Reading comprehension", "Backward digit recall")
)

cmi_cor_results <- do.call(
  rbind,
  lapply(cmi_cor_tests, function(v) {
    get_pearson_stats(
      data = beh_cmi,
      cohort_name = "CMI-HBN",
      xvar = v[[1]],
      yvar = v[[2]],
      x_label = v[[3]],
      y_label = v[[4]]
    )
  })
)

stanford_cor_results <- do.call(
  rbind,
  lapply(stanford_cor_tests, function(v) {
    get_pearson_stats(
      data = beh_su,
      cohort_name = "Stanford",
      xvar = v[[1]],
      yvar = v[[2]],
      x_label = v[[3]],
      y_label = v[[4]]
    )
  })
)

behavioral_cor_results <- rbind(cmi_cor_results, stanford_cor_results)

write.csv(
  behavioral_cor_results,
  "results/behavior/behavioral_correlation_results_full_statistics.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------
# 6. Group differences between CMI-HBN and Stanford
# ------------------------------------------------------------------

group_vars <- list(
  Age = list(x = beh_cmi$age, y = beh_su$age),
  FSIQ = list(x = beh_cmi$WISC_FSIQ, y = beh_su$fsiq),
  Numerical_operations = list(x = beh_cmi$numop_std, y = beh_su$numop_std),
  Math_reasoning = list(x = beh_cmi$mathprob_std, y = beh_su$mathrea_std),
  Word_reading = list(x = beh_cmi$wordread_std, y = beh_su$wordread_std),
  Reading_comprehension = list(x = beh_cmi$readcomp_std, y = beh_su$readcomp_std)
)

cont_group_results <- do.call(
  rbind,
  lapply(names(group_vars), function(v) {
    out <- get_cont_group_stats(
      x = group_vars[[v]]$x,
      y = group_vars[[v]]$y
    )
    out$Variable <- v
    out
  })
)

cont_group_results <- cont_group_results[, c(
  "Variable", "CMI", "Stanford", "Group_difference"
)]

# Sex distribution group comparison
gender_table <- data.frame(
  Dataset = c(rep("CMI-HBN", nrow(beh_cmi)), rep("Stanford", nrow(beh_su))),
  GenderF = c(beh_cmi$genderF, beh_su$genderF)
)

gender_table <- gender_table[complete.cases(gender_table), ]

gender_counts <- table(gender_table$Dataset, gender_table$GenderF)

chi_test <- chisq.test(gender_counts, correct = FALSE)
v <- cramers_v(gender_counts, ci = 0.95)

sex_group_result <- data.frame(
  Variable = "Sex, M/F",
  CMI = paste0(gender_counts["CMI-HBN", "0"], "/", gender_counts["CMI-HBN", "1"]),
  Stanford = paste0(gender_counts["Stanford", "0"], "/", gender_counts["Stanford", "1"]),
  Group_difference = sprintf(
    "\u03C7\u00B2(%d) = %.2f, p = %s, Cramer's V = %.2f, 95%% CI %s",
    round(unname(chi_test$parameter)),
    unname(chi_test$statistic),
    format_p(chi_test$p.value),
    v$Cramers_v,
    format_ci(v$CI_low, v$CI_high)
  ),
  stringsAsFactors = FALSE
)

behavioral_group_difference_results <- rbind(
  cont_group_results[cont_group_results$Variable == "Age", ],
  sex_group_result,
  cont_group_results[cont_group_results$Variable != "Age", ]
)

write.csv(
  behavioral_group_difference_results,
  "results/behavior/behavioral_group_difference_results_full_statistics.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------
# 7. Print results
# ------------------------------------------------------------------

cat("\nDescriptive statistics:\n")
print(behavioral_descriptives)

cat("\nSex counts:\n")
print(sex_counts)

cat("\nBehavioral correlation results:\n")
print(behavioral_cor_results)

cat("\nBehavioral group difference results:\n")
print(behavioral_group_difference_results)
