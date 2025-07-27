############################################################
# Identify Common Top 20% Brain Regions Across Cohorts
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script:
#   1. Loads Brainnetome (BN) atlas and Shirer network mappings.
#   2. Loads the top 20% ROI indices (based on CCA weights) 
#      for both math- and reading-related modes across CMI and Stanford.
#   3. Identifies ROIs common between the two cohorts.
#   4. Saves results as CSV files for both math and reading.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list = ls())

# Load required packages
library(readxl)

# ------------------------------------------------------------------
# 2. Load Brainnetome Atlas
# ------------------------------------------------------------------

bn_atlas_file <- "data/atlas/bn_atlas.xlsx"
bn_atlas <- read_excel(bn_atlas_file)

# Keep the first 218 regions (corresponding to cortical/subcortical ROIs)
bn_atlas <- bn_atlas[1:218, ]
bn_atlas$`Region ID` <- as.numeric(bn_atlas$`Region ID`)  # Ensure Region ID is numeric

# ------------------------------------------------------------------
# 3. Load Shirer-to-BN Mapping
# ------------------------------------------------------------------

bn_shier <- read.csv("data/atlas/shirer_to_bn.csv", header = FALSE)
colnames(bn_shier) <- c("BN_ID", "Shirer_ID", "Shirer_network")

# Keep only the first 218 regions
bn_shier <- bn_shier[1:218, ]
bn_shier$Shirer_network <- trimws(bn_shier$Shirer_network, which = "left")  # Remove leading spaces

# ------------------------------------------------------------------
# 4. Combine BN Atlas with Shirer Mapping
# ------------------------------------------------------------------

bn_network <- merge(bn_atlas, bn_shier, 
                    by.x = "Region ID", 
                    by.y = "BN_ID")

# Display the count of ROIs per Shirer network
table(bn_network$Shirer_network)

# ------------------------------------------------------------------
# 5. Load Top 20% ROI Indices for Math and Reading
# ------------------------------------------------------------------
# These CSV files contain ROI indices that correspond to the top 20% 
# of absolute CCA weights for each cohort and mode.

cmi_math <- read.csv("results/cca/brainmaps/xloading_m2_thr80_cmi_math_ageinmodel_SignChanged_cca.csv", header = FALSE)
su_math  <- read.csv("results/cca/brainmaps/xloading_m2_thr80_stanford_math_ageinmodel_SignChanged_cca.csv", header = FALSE)
cmi_read <- read.csv("results/cca/brainmaps/xloading_m2_thr80_cmi_reading_ageinmodel_SignChanged_cca.csv", header = FALSE)
su_read  <- read.csv("results/cca/brainmaps/xloading_m2_thr80_stanford_reading_ageinmodel_SignChanged_cca.csv", header = FALSE)

# ------------------------------------------------------------------
# 6. Identify Common Top ROIs Across Cohorts
# ------------------------------------------------------------------
# Find the intersection between CMI and Stanford top 20% ROIs for each mode.

df_math_comm <- subset(bn_network, `Region ID` %in% intersect(cmi_math$V1, su_math$V1))
df_read_comm <- subset(bn_network, `Region ID` %in% intersect(cmi_read$V1, su_read$V1))

# ------------------------------------------------------------------
# 7. Save Results
# ------------------------------------------------------------------

write.table(df_math_comm, "results/cca/brainmaps/top20_math_common.csv", 
            sep = ",", row.names = FALSE)
write.table(df_read_comm, "results/cca/brainmaps/top20_read_common.csv", 
            sep = ",", row.names = FALSE)

cat("Common top 20% ROI lists saved successfully!\n")
