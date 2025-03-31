setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(ggplot2)
library(readxl)

# provide path to files
bn_atlas_file = "data/bn_atlas.xlsx"
bn_atlas = read_excel(bn_atlas_file)
bn_atlas = bn_atlas[c(1:218),c(1:3)]
# get fine-grained gyrus
sp = strsplit(bn_atlas$Description,"_")
bn_atlas$subgyrus = ""
for(i in c(1:218)){
  bn_atlas$subgyrus[i] = sp[[i]][1]
}

# bn mapping to Shier
bn_shier = read.csv("scripts/shier_to_bn/shier_to_bn.csv", header = F)
colnames(bn_shier) = c("BN_ID","Shier_ID","Shier_network")
bn_shier = bn_shier[c(1:218),]


# neurotransmitter maps
fname_receptors = 'data/Neurotransmitter/receptor_data_bn246.csv'
receptors = read.csv(fname_receptors, header = F)

receptor_names = c('5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'A4B2', 'CB1',
                   'D1', 'D2', 'DAT', 'GABAa', 'H3', 'M1', 'mGluR5', 'MOR', 'NET',
                   'NMDA', 'VAChT')

# brain maps (cmi)
fname_cmi = 'results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef_cmi_n760.csv'
# fname_cmi = 'results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef_cmi_n760.csv'

# brain maps (stanford)
fname_su = 'results/cca/combined_CCA_PCA_roi_smri_brainnetome_math_age_coefs_stanford_n231.csv'
# fname_su = 'results/cca/combined_CCA_PCA_roi_smri_brainnetome_reading_age_coefs_stanford_n231.csv'

# load data
bmaps_cmi = read.csv(fname_cmi)
bmaps_su = read.csv(fname_su)

# organize data (change sign of brain maps)
df_cmi = data.frame(-1*bmaps_cmi[,4], receptors[1:218,])
df_su = data.frame(-1*bmaps_su[,4], receptors[1:218,])
colnames(df_cmi) = c('bmap', receptor_names)
colnames(df_cmi) = make.names(colnames(df_cmi))
colnames(df_su) = colnames(df_cmi)

df_cmi$gyrus = bn_atlas$subgyrus
df_cmi$cortex = bn_atlas$Gyrus
df_cmi$shier = bn_shier$Shier_network
df_su$gyrus = bn_atlas$subgyrus
df_su$cortex = bn_atlas$Gyrus
df_su$shier = bn_shier$Shier_network

# output
bmap_name = 'math_gmv'  #'reading_gmv' # 

i = 18
# for(i in c(1:19)){
  receptor = receptor_names[i]
  outputf = paste('figures_updated/Fig6_scatter_plot_directionality/cmi/scatter_',bmap_name,'_',receptor,'_NMDA_color_shier14.eps', sep="")
  
  ggplot(df_cmi, aes(x=bmap, y=df_cmi[,i+1])) + 
    geom_point(aes(color = shier), size = 8, alpha = 0.5) +
    geom_text(aes(label = shier), vjust = 0.4, hjust = 0.5, size = 3, check_overlap = TRUE) +
    ylab(receptor) +
    geom_smooth(method=lm, se=TRUE, color = "black") +
    theme_classic()  +
    theme(legend.position = "none")
  
  ggsave(outputf,units="in", width=8, height=8, device=cairo_ps)
  
  outputf = paste('figures_updated/Fig6_scatter_plot_directionality/su/scatter_',bmap_name,'_',receptor,'_NMDA_color_shier14.eps', sep="")
  
  ggplot(df_su, aes(x=bmap, y=df_su[,i+1])) + 
    geom_point(aes(color = shier), size = 8, alpha = 0.5) +
    geom_text(aes(label = shier), vjust = 0.4, hjust = 0.5, size = 3, check_overlap = TRUE) +
    ylab(receptor) +
    geom_smooth(method=lm, se=TRUE, color = "black") +
    theme_classic() +
    theme(legend.position = "none")
  
  ggsave(outputf,units="in", width=8, height=8, device=cairo_ps)
  
# }


 


