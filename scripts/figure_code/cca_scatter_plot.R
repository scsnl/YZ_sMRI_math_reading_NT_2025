setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(ggplot2)

# fpath = 'results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_CMI_full_math_n760/'
# fname = 'CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData'
# outputf = 'figures_updated/Fig2_combine_CCA_mathread_SignChanged/scatter_mathmode_cmi_n760.eps'

fpath = 'results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_CMI_full_reading_n760/'
fname = 'CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData'
outputf = 'figures_updated/Fig2_combine_CCA_mathread_SignChanged/scatter_readmode_cmi_n760.eps'

# fpath = 'results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_SU_math/'
# fname = 'CCA_PCA_roi_gmv_brainnetome_baseline_wiat2_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData'
# outputf = 'figures_updated/Fig2_combine_CCA_mathread_SignChanged/scatter_mathmode_stanford_n231.eps'

# fpath = 'results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_SU_reading/'
# fname = 'CCA_PCA_roi_gmv_brainnetome_baseline_wiat2_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData'
# outputf = 'figures_updated/Fig2_combine_CCA_mathread_SignChanged/scatter_readmode_stanford_n231.eps'

load(paste(fpath,fname,sep=""))
df = data.frame(brain = res$Cx[,2], beh = res$Cy[,2])

ggplot(df, aes(x=brain, y=beh)) + 
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method=lm, se=TRUE) +
  theme_classic()

ggsave(outputf,units="in", width=5, height=5, device=cairo_ps)
