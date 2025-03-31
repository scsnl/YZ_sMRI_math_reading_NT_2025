setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

# load libraries
# library(lmerTest)
# library(bestNormalize)
library(BayesFactor)

# neurotransmitter maps
fname_receptors = 'data/Neurotransmitter/receptor_data_bn246.csv'
receptors = read.csv(fname_receptors, header = F)

receptor_names = c('5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
                   'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
                   'NMDA', 'VAChT')

# brain maps (cmi)
# fname_cmi = 'results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef_cmi_n760.csv'
fname_cmi = 'results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef_cmi_n760.csv'

# brain maps (stanford)
# fname_su = 'results/cca/combined_CCA_PCA_roi_smri_brainnetome_math_age_coefs_stanford_n231.csv'
fname_su = 'results/cca/combined_CCA_PCA_roi_smri_brainnetome_reading_age_coefs_stanford_n231.csv'

# load data
bmaps_cmi = read.csv(fname_cmi)
bmaps_su = read.csv(fname_su)

# organize data (change sign of brain maps)
df_cmi = data.frame(-1*bmaps_cmi[,4], receptors[1:218,])
df_su = data.frame(-1*bmaps_su[,4], receptors[1:218,])
colnames(df_cmi) = c('bmap', receptor_names)
colnames(df_cmi) <- make.names(colnames(df_cmi))
colnames(df_su) = colnames(df_cmi)

# output
# bmap_name = 'math_gmv' 
bmap_name = 'reading_gmv'
outputf = paste('results/BF10/replication_bf10_',bmap_name,'_receptors.csv', sep="")


# BF
bf10_ratio1 = c()
bf10_ratio2 = c()
for (i in 2:ncol(df_cmi)) {
  formula <- as.formula(paste("bmap ~ ", colnames(df_cmi)[i], sep=""))
  # BF for original dataset
  BF1 = lmBF(formula, data = df_cmi)
  # BF for replication dataset
  BF2 = lmBF(formula, data = df_su)
  # BF total
  BF_total = lmBF(formula, data=rbind(df_cmi,df_su))
  
  # ratio1: use cmi as original and su as repliaction
  ratio1 = extractBF(BF_total, logbf = F, onlybf = T) / extractBF(BF1, logbf = F, onlybf = T)
  # ratio2: use su as original and cmi as replication
  ratio2 = extractBF(BF_total, logbf = F, onlybf = T) / extractBF(BF2, logbf = F, onlybf = T)
  
  bf10_ratio1 = c(bf10_ratio1, ratio1)
  bf10_ratio2 = c(bf10_ratio2, ratio2)
}

res = data.frame(receptor = receptor_names, BF10_cmi_su = bf10_ratio1, BF10_su_cmi = bf10_ratio2)
res$bmap = bmap_name
idx = c(9:11,15,18,12,1:7,14,19,17,8,13,16)
res = res[idx,c(4,1,2,3)]
 
write.csv(res, file = outputf, row.names=F)

