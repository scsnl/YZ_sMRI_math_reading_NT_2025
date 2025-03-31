setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(readxl)

# cmi cohort
beh_cmi = read.csv('data/subjectlists/cmi/2024_04_22_CMI_supercleanTD_completeNP_7-13y_n858_withNPinNeed_updatename.csv')
gmv_file = read.csv('data/gmv_lCT_cmi_n760.csv')
beh_cmi = subset(beh_cmi, oak_id %in% gmv_file$pid)

cor.test(beh_cmi$numop_std, beh_cmi$wordread_std)
cor.test(beh_cmi$numop_std, beh_cmi$readcomp_std)
cor.test(beh_cmi$mathprob_std, beh_cmi$wordread_std)
cor.test(beh_cmi$mathprob_std, beh_cmi$readcomp_std)
cor.test(beh_cmi$age, beh_cmi$numop_std)
cor.test(beh_cmi$age, beh_cmi$mathprob_std)
cor.test(beh_cmi$age, beh_cmi$wordread_std)
cor.test(beh_cmi$age, beh_cmi$readcomp_std)

# stanford cohort
beh_su = read.csv("data/subjectlists/stanford/subjectlist_baseline_wiat2_covariates_woEO.csv")
idx = which(beh_su$pid %in% c(203, 209, 7401, 7698))
beh_su = beh_su[-idx,]

cor.test(beh_su$numop_std, beh_su$wordread_std)
cor.test(beh_su$numop_std, beh_su$readcomp_std)
cor.test(beh_su$mathrea_std, beh_su$wordread_std)
cor.test(beh_su$mathrea_std, beh_su$readcomp_std)
cor.test(beh_su$age, beh_su$numop_std)
cor.test(beh_su$age, beh_su$mathrea_std)
cor.test(beh_su$age, beh_su$wordread_std)
cor.test(beh_su$age, beh_su$readcomp_std)

# group differences
t.test(beh_cmi$age, beh_su$age, var.equal=T)
t.test(beh_cmi$WISC_FSIQ, beh_su$fsiq, var.equal=T)
t.test(beh_cmi$numop_std, beh_su$numop_std, var.equal=T)
t.test(beh_cmi$mathprob_std, beh_su$mathrea_std, var.equal=T)
t.test(beh_cmi$wordread_std, beh_su$wordread_std, var.equal=T)
t.test(beh_cmi$readcomp_std, beh_su$readcomp_std, var.equal=T)

beh_cmi$genderF = 1
beh_cmi$genderF[which(beh_cmi$Sex == 'M')] = 0
# Create a contingency table of genderF for each dataset
gender_table <- data.frame(
  Dataset = c(rep("CMI", nrow(beh_cmi)), rep("SU", nrow(beh_su))),
  GenderF = c(beh_cmi$genderF, beh_su$genderF)
)

# Summarize the counts for each gender in both datasets
gender_counts <- table(gender_table$Dataset, gender_table$GenderF)
# Perform chi-squared test to see if there's a significant difference in sex between the two datasets
chi_test <- chisq.test(gender_counts)
chi_test
