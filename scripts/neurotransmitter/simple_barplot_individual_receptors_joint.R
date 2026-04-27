setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")
rm(list=ls())

library(ggplot2) 
library(forcats)

# CMI results (shared mode)
fname = "results/neurotransmitter/cmi/shared/individual_neurotransmitter_regression_results_shared_mode2_wFDRp.csv"
postfix = "barplot_adj_r2_shared_cmi"
output_path = 'results/neurotransmitter/cmi/shared'
outputf = paste(output_path, '/', postfix, '.eps',sep="")

# # Stanford results (shared mode)
# fname = "results/neurotransmitter/stanford/shared/individual_neurotransmitter_regression_results_shared_mode2_wFDRp.csv"
# postfix = "barplot_adj_r2_shared_stanford"
# output_path = 'results/neurotransmitter/stanford/shared'
# outputf = paste(output_path, '/', postfix, '.eps',sep="")

df = read.csv(fname)
df$receptor_set = as.factor(df$receptor_set)
df$receptor_set = factor(df$receptor_set, 
                         levels=levels(df$receptor_set)[c(9:11,15,18,12,1:6,7,14,19,17,8,13,16)])

ggplot(df, aes(x=fct_rev(receptor_set), y=adj_r2)) + 
  geom_bar(stat = "identity",width=0.86,fill="darkorange") +
  ylim(c(-0.01,0.6)) +
  coord_flip() +
  theme_classic() +
  geom_text(aes(label=sprintf("p = %.3f", fdrp)),
            hjust=-0.1, # Adjusts the horizontal position of the text
            vjust=0.5, # Adjusts the vertical position of the text
            color="black", # Text color
            size=5) +
  theme(axis.text.x = element_text(face = "bold", size = 14),  # Bold and increase font size of x-axis labels
        axis.text.y = element_text(face = "bold", size = 14),  # Bold and increase font size of y-axis labels
        axis.title.x = element_text(size = 14, face = "bold"),  # Increase font size of x-axis title
        axis.title.y = element_text(size = 14, face = "bold")) 
ggsave(outputf,units="in", width=5, height=4)

