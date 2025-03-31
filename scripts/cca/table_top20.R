setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(ggplot2)
library(readxl)

# provide path to files
bn_atlas_file = "data/bn_atlas.xlsx"
bn_atlas = read_excel(bn_atlas_file)
bn_atlas = bn_atlas[c(1:218),]
bn_atlas$`Region ID` = as.numeric(bn_atlas$`Region ID`)

# bn mapping to Shier
bn_shier = read.csv("scripts/shier_to_bn/shier_to_bn.csv", header = F)
colnames(bn_shier) = c("BN_ID","Shier_ID","Shier_network")
bn_shier = bn_shier[c(1:218),]
bn_shier$Shier_network = trimws(bn_shier$Shier_network, which = "left")

# combind
bn_network = merge(bn_atlas, bn_shier, by.x = "Region ID", by.y = "BN_ID")
table(bn_network$Shier_network)

# load top 20% regions
cmi_math = read.csv("figures_updated/FigS1_top_regions/top20/xloading_m2_math_gmv_thr80_SignChanged_cca_cmi_n760.csv",header=F)
su_math = read.csv("figures_updated/FigS1_top_regions/top20/xloading_m2_math_gmv_thr80_SignChanged_cca_stanford_n231.csv",header=F)
cmi_read = read.csv("figures_updated/FigS1_top_regions/top20/xloading_m2_read_gmv_thr80_SignChanged_cca_cmi_n760.csv",header=F)
su_read = read.csv("figures_updated/FigS1_top_regions/top20/xloading_m2_read_gmv_thr80_SignChanged_cca_stanford_n231.csv",header=F)

# merge
df_cmi_math = subset(bn_network, `Region ID` %in% cmi_math$V1)
df_su_math = subset(bn_network, `Region ID` %in% su_math$V1)
df_math_comm = subset(bn_network, `Region ID` %in% intersect(cmi_math$V1, su_math$V1))

df_cmi_read = subset(bn_network, `Region ID` %in% cmi_read$V1)
df_su_read = subset(bn_network, `Region ID` %in% su_read$V1)
df_read_comm = subset(bn_network, `Region ID` %in% intersect(cmi_read$V1, su_read$V1))

write.table(df_math_comm, "figures_updated/FigS1_top_regions/top20/top20_math_common.csv", sep=",", row.names=FALSE)
write.table(df_read_comm, "figures_updated/FigS1_top_regions/top20/top20_read_common.csv", sep=",", row.names=FALSE)

