setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

# ------------------------------------------------------------
# Get R package versions for reporting summary
# ------------------------------------------------------------

packages <- c(
  "CCA",
  "CCP",
  "ppcor",
  "BayesFactor",
  "permute",
  "psych",
  "readxl",
  "R.matlab",
  "reshape2",
  "ggplot2",
  "ggrepel",
  "ggseg"
)

get_pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(packageVersion(pkg))
  } else {
    NA_character_
  }
}

r_versions <- data.frame(
  software_or_package = c("R", packages),
  version = c(
    as.character(getRversion()),
    sapply(packages, get_pkg_version)
  ),
  stringsAsFactors = FALSE
)

print(r_versions)

write.csv(
  r_versions,
  "R_package_versions.csv",
  row.names = FALSE
)

# Optional: save full session information
sink("R_sessionInfo.txt")
sessionInfo()
sink()