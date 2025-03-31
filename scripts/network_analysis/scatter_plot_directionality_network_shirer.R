setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(ggplot2)
library(readxl)

# provide path to files
bn_atlas_file = "data/bn_atlas.xlsx"
bn_atlas = read_excel(bn_atlas_file)
bn_atlas = bn_atlas[c(1:218),c(1:3)]
bn_atlas$`Region ID` = as.numeric(bn_atlas$`Region ID`)
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
bn_shier$Shier_network = trimws(bn_shier$Shier_network, which = "left")


# combind
bn_network = merge(bn_atlas, bn_shier, by.x = "Region ID", by.y = "BN_ID")
table(bn_network$Shier_network)

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

# output
bmap_name = 'math_gmv'  
# bmap_name = 'reading_gmv' 


################################
# load data
bmaps_cmi = read.csv(fname_cmi)
bmaps_su = read.csv(fname_su)

# organize data (change sign of brain maps)
df_cmi = data.frame(-1*bmaps_cmi[,4], receptors[1:218,])
df_su = data.frame(-1*bmaps_su[,4], receptors[1:218,])
colnames(df_cmi) = c('bmap', receptor_names)
colnames(df_cmi) <- make.names(colnames(df_cmi))
colnames(df_su) = colnames(df_cmi)

###################
df = df_cmi
cohort = "cmi"
# df = df_su
# cohort = "su"

df = df[,c(1,which(receptor_names == "NMDA")+1)]
# df[,1] = scale(df[,1])
# df[,2] = scale(df[,2])
df$gyrus = bn_network$Gyrus
df$subgyrus = bn_network$subgyrus
df$shier = bn_network$Shier_network

###################
partition = c("gyrus","subgyrus","shier")
i = 3
# for(i in c(1:4)){
  # output_csv = paste("results/NMDA_directionality_scale/nmda_", bmap_name, "_", cohort, "_", partition[i], ".csv",sep="")
  output_csv = paste("results/NMDA_directionality/nmda_", bmap_name, "_", cohort, "_", partition[i], ".csv",sep="")

  networks = unique(df[,2+i])
  nn = length(networks) # number of networks
  
  res = data.frame(network=networks, t=numeric(nn), df=numeric(nn), r=numeric(nn), p=numeric(nn))
  # for each network, get the correlation
  for (j in c(1:nn)){
    df_net = subset(df, df[,2+i] == networks[j])
    
    # stats
    tmp = cor.test(df_net[,1], df_net[,2])
    res$t[j] = tmp$statistic
    res$df[j] = tmp$parameter
    res$r[j] = tmp$estimate
    res$p[j] = tmp$p.value
    
    # plots
    ylabel = paste(colnames(df_net)[2]," (", networks[j], ")",sep="")
    xlabel = paste("brain maps of ", bmap_name, " (", networks[j], ")",sep="")
      
    ggplot(df_net, aes(x=bmap, y=df_net[,2])) +
      geom_point(size = 5, alpha = 0.5) +
      ylab(ylabel) +
      xlab(xlabel) +
      geom_smooth(method=lm, se=TRUE) +
      theme_classic()  +
      theme(legend.position = "none")
    
    # output_fig = paste("results/NMDA_directionality_scale/nmda_", bmap_name, "_", cohort, "_", partition[i], "_", networks[j],".eps",sep="")
    output_fig = paste("results/NMDA_directionality/nmda_", bmap_name, "_", cohort, "_", partition[i], "_", networks[j],".eps",sep="")
    ggsave(output_fig, units="in", width=5, height=5, device=cairo_ps)
  }
  res$fdrp = p.adjust(res$p, method="fdr", n=length(res$p))
  write.csv(res, file=output_csv, row.names = F)
# }
