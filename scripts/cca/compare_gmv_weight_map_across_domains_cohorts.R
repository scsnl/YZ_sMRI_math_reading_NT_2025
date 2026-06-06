############################################################
# Compare GMV weight maps across cohorts/domains
# Rigorous comparison using cocor + figures
# Author: Yuan Zhang
# Date: 2026-03-24
############################################################


rm(list = ls())

library(Hmisc)
library(corrplot)
library(ggplot2)
library(cocor)

# ----------------------------------------------------------
# 1. File paths
# ----------------------------------------------------------
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

cmi_math_file  <- "results/sex_analysis/cmi/gmv_weights_sex_math.csv"
cmi_read_file  <- "results/sex_analysis/cmi/gmv_weights_sex_read.csv"
stan_math_file <- "results/sex_analysis/stanford/gmv_weights_sex_math.csv"
stan_read_file <- "results/sex_analysis/stanford/gmv_weights_sex_read.csv"

out_dir <- "results/gmv_map_crossdomain_check/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 2. Read maps
# Assumption: full-sample map is in column "Full"
# ----------------------------------------------------------
cmi_math  <- read.csv(cmi_math_file)
cmi_read  <- read.csv(cmi_read_file)
stan_math <- read.csv(stan_math_file)
stan_read <- read.csv(stan_read_file)

map_df <- data.frame(
  CMI_math       = cmi_math$Full,
  CMI_reading    = cmi_read$Full,
  Stanford_math  = stan_math$Full,
  Stanford_read  = stan_read$Full
)

# sanity check
stopifnot(nrow(map_df) > 3)

# ----------------------------------------------------------
# 3. Pairwise correlations
# ----------------------------------------------------------
rc <- rcorr(as.matrix(map_df), type = "pearson")

cor_mat <- rc$r
p_mat   <- rc$P

write.csv(cor_mat, file.path(out_dir, "gmv_map_correlation_matrix.csv"))
write.csv(p_mat,   file.path(out_dir, "gmv_map_pvalue_matrix.csv"))

pairs <- combn(colnames(map_df), 2, simplify = FALSE)
pairwise_results <- do.call(rbind, lapply(pairs, function(pair) {
  x <- map_df[[pair[1]]]
  y <- map_df[[pair[2]]]
  idx <- which(!is.na(x) & !is.na(y))
  test <- cor.test(x[idx], y[idx], method = "pearson")
  
  data.frame(
    var1 = pair[1],
    var2 = pair[2],
    n = length(idx),
    r = unname(test$estimate),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    ci_lower = test$conf.int[1],
    ci_upper = test$conf.int[2],
    stringsAsFactors = FALSE
  )
}))
write.csv(pairwise_results,
          file.path(out_dir, "gmv_map_pairwise_correlations.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 4. Rigorous comparison of dependent overlapping correlations
# ----------------------------------------------------------
# Comparison 1:
#   Is cor(CMI_math, Stanford_math) > cor(CMI_math, CMI_reading)?
r_jk_1 <- cor(map_df$CMI_math, map_df$Stanford_math)   # CMI_math vs Stanford_math
r_jh_1 <- cor(map_df$CMI_math, map_df$CMI_reading)     # CMI_math vs CMI_reading
r_kh_1 <- cor(map_df$Stanford_math, map_df$CMI_reading)
n1 <- nrow(map_df)

comp1 <- cocor.dep.groups.overlap(
  r.jk = r_jk_1,
  r.jh = r_jh_1,
  r.kh = r_kh_1,
  n = n1,
  alternative = "two.sided",
  test = c("pearson1898", "hotelling1940", "williams1959", "steiger1980", "zou2007")
)

# Comparison 2:
#   Is cor(CMI_reading, Stanford_read) > cor(CMI_reading, CMI_math)?
r_jk_2 <- cor(map_df$CMI_reading, map_df$Stanford_read)  # CMI_reading vs Stanford_read
r_jh_2 <- cor(map_df$CMI_reading, map_df$CMI_math)       # CMI_reading vs CMI_math
r_kh_2 <- cor(map_df$Stanford_read, map_df$CMI_math)
n2 <- nrow(map_df)

comp2 <- cocor.dep.groups.overlap(
  r.jk = r_jk_2,
  r.jh = r_jh_2,
  r.kh = r_kh_2,
  n = n2,
  alternative = "two.sided",
  test = c("pearson1898", "hotelling1940", "williams1959", "steiger1980", "zou2007")
)

# Save cocor text output
sink(file.path(out_dir, "gmv_map_cocor_results.txt"))
cat("====================================================\n")
cat("Comparison 1: CMI_math-Stanford_math vs CMI_math-CMI_reading\n")
cat("====================================================\n")
print(comp1)
cat("\n\n")
cat("====================================================\n")
cat("Comparison 2: CMI_reading-Stanford_read vs CMI_reading-CMI_math\n")
cat("====================================================\n")
print(comp2)
sink()

# ----------------------------------------------------------
# 5. Heatmap figure
# ----------------------------------------------------------
plot_corr <- function(cor_matrix, outfile) {
  my_col <- colorRampPalette(c("white", "orange", "red4"))(200)
  
  postscript(outfile, width = 5.5, height = 5.5, horizontal = FALSE, onefile = FALSE, paper = "special")
  corrplot(
    cor_matrix,
    method = "color",
    col = my_col,
    type = "upper",
    diag = FALSE,
    addCoef.col = "black",
    number.cex = 0.8,
    tl.col = "black",
    tl.srt = 35,
    cl.pos = "r",
    is.corr = TRUE,
    mar = c(0, 0, 1, 0)
  )
  dev.off()
}

plot_corr(cor_mat, file.path(out_dir, "gmv_map_correlation_heatmap.eps"))

# ----------------------------------------------------------
# 6. Scatter plots for the two key comparisons
# ----------------------------------------------------------
make_scatter <- function(df, xvar, yvar, title_text, outfile) {
  test <- cor.test(df[[xvar]], df[[yvar]], method = "pearson")
  
  p <- ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(size = 2.2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "tomato") +
    annotate(
      "text",
      x = min(df[[xvar]], na.rm = TRUE),
      y = max(df[[yvar]], na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      size = 4,
      label = paste0(
        "r = ", sprintf("%.3f", unname(test$estimate)),
        "\np = ", formatC(test$p.value, format = "e", digits = 2)
      )
    ) +
    labs(
      title = title_text,
      x = xvar,
      y = yvar
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(outfile, p, width = 5.2, height = 4.8)
}

make_scatter(
  map_df,
  "CMI_math", "Stanford_math",
  "CMI math vs Stanford math GMV maps",
  file.path(out_dir, "scatter_CMI_math_vs_Stanford_math.png")
)

make_scatter(
  map_df,
  "CMI_math", "CMI_reading",
  "CMI math vs CMI reading GMV maps",
  file.path(out_dir, "scatter_CMI_math_vs_CMI_reading.png")
)

make_scatter(
  map_df,
  "CMI_reading", "Stanford_read",
  "CMI reading vs Stanford reading GMV maps",
  file.path(out_dir, "scatter_CMI_reading_vs_Stanford_read.png")
)

make_scatter(
  map_df,
  "CMI_reading", "CMI_math",
  "CMI reading vs CMI math GMV maps",
  file.path(out_dir, "scatter_CMI_reading_vs_CMI_math.png")
)

# ----------------------------------------------------------
# 7. Simple summary table
# ----------------------------------------------------------
get_cor_stats <- function(df, xvar, yvar, comparison_name) {
  x <- df[[xvar]]
  y <- df[[yvar]]
  idx <- which(!is.na(x) & !is.na(y))
  test <- cor.test(x[idx], y[idx], method = "pearson")
  
  data.frame(
    comparison = comparison_name,
    var1 = xvar,
    var2 = yvar,
    n = length(idx),
    r = unname(test$estimate),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    ci_lower = test$conf.int[1],
    ci_upper = test$conf.int[2],
    stringsAsFactors = FALSE
  )
}

summary_df <- rbind(
  get_cor_stats(
    map_df,
    "CMI_math", "Stanford_math",
    "CMI_math_vs_Stanford_math"
  ),
  get_cor_stats(
    map_df,
    "CMI_math", "CMI_reading",
    "CMI_math_vs_CMI_reading"
  ),
  get_cor_stats(
    map_df,
    "CMI_reading", "Stanford_read",
    "CMI_reading_vs_Stanford_read"
  ),
  get_cor_stats(
    map_df,
    "CMI_reading", "CMI_math",
    "CMI_reading_vs_CMI_math"
  )
)

write.csv(summary_df,
          file.path(out_dir, "gmv_map_key_correlation_summary.csv"),
          row.names = FALSE)

# ----------------------------------------------------------
# 8. Console output
# ----------------------------------------------------------
cat("\n==============================\n")
cat("Key correlations\n")
cat("==============================\n")
print(summary_df)

cat("\nCondition A:\n")
cat("CMI_math vs Stanford_math > CMI_math vs CMI_reading ? ",
    r_jk_1 > r_jh_1, "\n")

cat("Condition B:\n")
cat("CMI_reading vs Stanford_read > CMI_reading vs CMI_math ? ",
    r_jk_2 > r_jh_2, "\n")

cat("\nResults saved to:\n", out_dir, "\n")

cat("\n==============================\n")
cat("Formatted correlation results\n")
cat("==============================\n")

for (i in 1:nrow(summary_df)) {
  cat(
    summary_df$comparison[i], ": ",
    "r = ", sprintf("%.3f", summary_df$r[i]),
    ", t(", summary_df$df[i], ") = ", sprintf("%.3f", summary_df$t[i]),
    ", p = ", formatC(summary_df$p[i], format = "e", digits = 2),
    ", 95% CI [", sprintf("%.3f", summary_df$ci_lower[i]), ", ",
    sprintf("%.3f", summary_df$ci_upper[i]), "]",
    ", n = ", summary_df$n[i],
    "\n",
    sep = ""
  )
}