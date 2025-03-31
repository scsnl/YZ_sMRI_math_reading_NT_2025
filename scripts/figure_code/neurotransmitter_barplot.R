setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(ggplot2) 
library(forcats)

# # cmi math
# fname = "results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_individual_neurotransmitter_regression_results_cmi_n760.csv"
# # cmi reading
# fname = "results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_individual_neurotransmitter_regression_results_cmi_n760.csv"
# # stanford math
# fname = "results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_individual_neurotransmitter_regression_results_stanford_n231.csv"
# stanford reading
fname = "results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_individual_neurotransmitter_regression_results_stanford_n231.csv"


# postfix = "barplot_adj_r2_math_gmv_individual_neurotransmitter_cmi_n760"
# postfix = "barplot_adj_r2_reading_gmv_individual_neurotransmitter_cmi_n760"
# postfix = "barplot_adj_r2_math_gmv_individual_neurotransmitter_stanford_n231"
postfix = "barplot_adj_r2_reading_gmv_individual_neurotransmitter_stanford_n231"

df = read.csv(fname)
rows = seq(from = 2, to = 38, by = 2)
df = df[rows,]

output_path = 'figures_updated/Fig4_individual_neurotransmitters_cmi'
outputf = paste(output_path, '/', postfix, '.eps',sep="")

df$receptor_set = as.factor(df$receptor_set)
df$receptor_set = factor(df$receptor_set, 
                         levels=levels(df$receptor_set)[c(9:11,15,18,12,1:6,7,14,19,17,8,13,16)])

ggplot(df, aes(x=fct_rev(receptor_set), y=adj_r2)) + 
  geom_bar(stat = "identity",width=0.86,fill="darkorange") +
  ylim(c(-0.01,0.6)) +
  coord_flip() +
  theme_classic() +
  # theme_bw() +
  # geom_text(aes(label=sprintf("p = %.3f", moran_p)), 
  #           hjust=-0.1, # Adjusts the horizontal position of the text
  #           vjust=0.5, # Adjusts the vertical position of the text
  #           color="black", # Text color
  #           size=5) +
  theme(axis.text.x = element_text(face = "bold", size = 14),  # Bold and increase font size of x-axis labels
        axis.text.y = element_text(face = "bold", size = 14),  # Bold and increase font size of y-axis labels
        axis.title.x = element_text(size = 14, face = "bold"),  # Increase font size of x-axis title
        axis.title.y = element_text(size = 14, face = "bold")) 
ggsave(outputf,units="in", width=5, height=4)

