setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(BayesFactor)

# neurotransmitter maps
fname_receptors = 'data/Neurotransmitter/receptor_data_bn246.csv'
receptors = read.csv(fname_receptors, header = F)

receptor_names = c('5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
                  'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
                  'NMDA', 'VAChT')

# brain maps (cmi)
fname = 'results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef_cmi_n760.csv'
# fname = 'results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef_cmi_n760.csv'

# brain maps (stanford)
# fname = 'results/cca/combined_CCA_PCA_roi_smri_brainnetome_math_age_coefs_stanford_n231.csv'
# fname = 'results/cca/combined_CCA_PCA_roi_smri_brainnetome_reading_age_coefs_stanford_n231.csv'


bmaps = read.csv(fname)

# organize data (change sign of brain maps)
df = data.frame(-1*bmaps[,4], receptors[1:218,])
colnames(df) = c('bmap', receptor_names)
colnames(df) <- make.names(colnames(df))

# output
bmap_name = 'math_gmv'
# bmap_name = 'reading_gmv'
outputf = paste('results/BF10/bf10_',bmap_name,'_receptors_cmi_n760.csv', sep="")
# outputf = paste('results/BF10/bf10_',bmap_name,'_receptors_stanford_n231.csv', sep="")

# BF
bf10 = c()
for (i in 2:ncol(df)) {
  formula <- as.formula(paste("bmap ~ ", colnames(df)[i], sep=""))
  bf = lmBF(formula, data = df)
  bf_value = exp(bf@bayesFactor$bf)
  bf10 = c(bf10, bf_value)
}

res = data.frame(receptor = receptor_names, BF10 = bf10)
res$bmap = bmap_name
idx = c(9:11,15,18,12,1:7,14,19,17,8,13,16)
res = res[idx,c(3,1,2)]

write.csv(res, file = outputf, row.names=F)



