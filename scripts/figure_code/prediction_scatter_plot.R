setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(ggplot2)
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
bmap_name = 'reading_gmv' # 'math_gmv' 
outputf = paste('results/BF10/replication_bf10_',bmap_name,'_receptors.csv', sep="")

ggplot(df, aes(x=predicted_brain_score, y=predicted_beh_score)) + 
  geom_point(size = 3, alpha = 0.5) +
  # geom_smooth(method=lm, se=FALSE) +
  geom_smooth(method=lm, se=TRUE) +
  theme_classic()

ggsave(outputf,units="in", width=5, height=5, device=cairo_ps)
