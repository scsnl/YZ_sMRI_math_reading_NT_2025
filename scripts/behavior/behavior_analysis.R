############################################################
# Behavioral Analysis Script
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script performs basic behavioral data analysis for 
# two cohorts (CMI and Stanford). It computes correlations 
# between behavioral measures, compares group differences 
# (t-tests), and examines gender distributions (chi-squared test).
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory (update this path to your local environment)
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear all objects from the environment
rm(list = ls())

# ------------------------------------------------------------------
# 2. Load Data
# ------------------------------------------------------------------

# Load behavioral data for CMI cohort (n = 760)
beh_cmi <- read.csv('data/subjectlist/subjectlist_cmi_n760.csv')

# ------------------------------------------------------------------
# 3. Correlation Analysis (CMI Cohort)
# ------------------------------------------------------------------
# Compute Pearson correlations between math, reading, and age.
# _std = grade_normed scores
cor.test(beh_cmi$numop_std, beh_cmi$wordread_std)   # Number operations vs. word reading
cor.test(beh_cmi$numop_std, beh_cmi$readcomp_std)   # Number operations vs. reading comprehension
cor.test(beh_cmi$mathprob_std, beh_cmi$wordread_std) # Math problem-solving vs. word reading
cor.test(beh_cmi$mathprob_std, beh_cmi$readcomp_std) # Math problem-solving vs. reading comprehension

# Correlations with age
cor.test(beh_cmi$age, beh_cmi$numop_std)
cor.test(beh_cmi$age, beh_cmi$mathprob_std)
cor.test(beh_cmi$age, beh_cmi$wordread_std)
cor.test(beh_cmi$age, beh_cmi$readcomp_std)

# ------------------------------------------------------------------
# 4. Correlation Analysis (Stanford Cohort)
# ------------------------------------------------------------------

# Load behavioral data for Stanford cohort (n = 231)
beh_su <- read.csv("data/subjectlist/subjectlist_stanford_n231.csv")

# Compute correlations between math, reading, and age
cor.test(beh_su$numop_std, beh_su$wordread_std)
cor.test(beh_su$numop_std, beh_su$readcomp_std)
cor.test(beh_su$mathrea_std, beh_su$wordread_std)   # mathrea_std = math reasoning
cor.test(beh_su$mathrea_std, beh_su$readcomp_std)

# Correlations with age
cor.test(beh_su$age, beh_su$numop_std)
cor.test(beh_su$age, beh_su$mathrea_std)
cor.test(beh_su$age, beh_su$wordread_std)
cor.test(beh_su$age, beh_su$readcomp_std)

# ------------------------------------------------------------------
# 5. Group Differences (CMI vs. Stanford)
# ------------------------------------------------------------------
# Independent-samples t-tests for age, IQ, and standardized measures.

t.test(beh_cmi$age, beh_su$age, var.equal = TRUE)
t.test(beh_cmi$WISC_FSIQ, beh_su$fsiq, var.equal = TRUE)
t.test(beh_cmi$numop_std, beh_su$numop_std, var.equal = TRUE)
t.test(beh_cmi$mathprob_std, beh_su$mathrea_std, var.equal = TRUE)
t.test(beh_cmi$wordread_std, beh_su$wordread_std, var.equal = TRUE)
t.test(beh_cmi$readcomp_std, beh_su$readcomp_std, var.equal = TRUE)

# ------------------------------------------------------------------
# 6. Gender Distribution Analysis
# ------------------------------------------------------------------

# Create a binary gender variable for the CMI cohort (1 = Female, 0 = Male)
beh_cmi$genderF <- 1
beh_cmi$genderF[which(beh_cmi$Sex == 'M')] <- 0

# Combine datasets for gender analysis
gender_table <- data.frame(
  Dataset = c(rep("CMI", nrow(beh_cmi)), rep("SU", nrow(beh_su))),
  GenderF = c(beh_cmi$genderF, beh_su$genderF)
)

# Summarize counts for each gender by dataset
gender_counts <- table(gender_table$Dataset, gender_table$GenderF)

# Perform chi-squared test for gender distribution differences
chi_test <- chisq.test(gender_counts)
chi_test
