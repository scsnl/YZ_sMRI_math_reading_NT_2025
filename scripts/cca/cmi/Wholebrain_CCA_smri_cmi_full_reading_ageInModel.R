setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())
source("scripts/newCCA_wholebrain_cmi/myRFunc_updated_smri_cmi_20240718.R")

#############################
# provide path to files
bn_atlas_file = "data/bn_atlas.xlsx"
beh_file = 'data/subjectlists/cmi/2024_04_22_CMI_supercleanTD_completeNP_7-13y_n858_withNPinNeed_updatename.csv'
gmv_file = 'data/gmv_lCT_cmi_n760.csv'

# apply PCA
apply_pca = 1
prefix = "CCA_PCA"
# apply_pca = 0
# prefix = "CCA"

# output path
output_path = "results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_CMI_full_reading_n760"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
condition = "roi_gmv_brainnetome"

postfix = "readstd_ageinmodel"
cols = c("age","wordread_std","readcomp_std")
use_resid = 0

########################
### data preparation ###
########################
# load data
data = load_data(bn_atlas_file, beh_file, gmv_file)
idx = which(data$beh$oak_id %in% data$pids)
data$beh = data$beh[idx,]
sum(data$beh$oak_id == data$pids) # checking brain and behavior mapping

# check data completeness
dim(data$brain)
sum(complete.cases(data$brain))
sum(complete.cases(data$beh))

# check extreme outliers for brain and behavioral data
# Function to identify outliers using IQR
identify_outliers <- function(column) {
  Q1 <- quantile(column, 0.25)
  Q3 <- quantile(column, 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 3 * IQR
  upper_bound <- Q3 + 3 * IQR
  which(column < lower_bound | column > upper_bound)
}
# identify extreme outliers across all brain regions
outliers <- lapply(as.data.frame(data$brain), identify_outliers)
brain_exo_idx = sort(unique(unlist(outliers))) 

# identify behavioral outliers
beh_df = as.data.frame(data$beh)
col_of_interest = c("age","mathprob_std","numop_std","readcomp_std","wordread_std")
beh_df = beh_df[,which(colnames(beh_df) %in% col_of_interest)]
outliers <- lapply(beh_df, identify_outliers)
beh_exo_idx = sort(unique(unlist(outliers)))

exo_idx = unique(c(beh_exo_idx, brain_exo_idx))
  
# remove extreme outliers
data$beh = data$beh[-exo_idx,]
data$brain = data$brain[-exo_idx,]
data$pids = data$pids[-exo_idx]
data$numsub = data$numsub - length(exo_idx)

# apply PCA
if(apply_pca == 1 & use_resid == 1){
  my.pca = mypca(data_resid$brain_resid, 0.9)
  Xorig = data_resid$brain_resid
}else if(apply_pca == 1 & use_resid == 0){
  my.pca = mypca(data$brain, 0.9)
  Xorig = data$brain
}else if(apply_pca == 0 & use_resid == 0){
  Xorig = data$brain
}else if(apply_pca == 0 & use_resid == 1){
  Xorig = data_resid$brain_resid
}


########################
######### CCA ##########
########################
if(apply_pca == 1){
  X = my.pca$x[,c(1:my.pca$numPC)]
}else{
  if(use_resid == 1){
    X = data_resid$brain_resid
  }else{
    X = data$brain 
  }
}

if(use_resid == 1){
  Y = data_resid$beh[,which(colnames(data_resid$beh) %in% cols)]
}else{
  Y = data$beh[,which(colnames(data$beh) %in% cols)]
}

# res includes canonical correlations
res = mycca_noscale(X,Y,5000)

# CCA model significance
res$Pillai
res$pillai.p
res$pillai.p.perm

# canonical correlations
res$cancor
res$cancor^2
# canonical correlations significance
res$cancor.p.perm
res$cancor.p.perm.adjust

# dimension test
p.asym(res$cancor, nrow(X), ncol(X), ncol(Y), tstat = "Wilks")

# save results
if(apply_pca == 1 & use_resid == 1){
  fname = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_perm5000_pcaNoscale_ccaNoScale.RData',sep="")
  save(data, data_resid, Xorig, X, Y, res, my.pca, file = fname)
}else if(apply_pca == 1 & use_resid == 0){
  fname = paste(output_path,'/', prefix, '_',condition, '_',postfix,'_perm5000_pcaNoscale_ccaNoScale.RData',sep="")
  save(data, Xorig, X, Y, res, my.pca, file = fname)
}else if(apply_pca == 0 & use_resid == 1){
  fname = paste(output_path,'/',prefix, '_',condition, '_',postfix,'_perm5000_ccaNoScale.RData',sep="")
  save(data, data_resid, Xorig, X, Y, res,  file = fname)
}else{
  fname = paste(output_path,'/',prefix, '_',condition, '_',postfix,'_perm5000_ccaNoScale.RData',sep="")
  save(data, Xorig, X, Y, res,  file = fname)
}


### 
# writting weights
write_weights(res, Xorig, Y, apply_pca, output_path, prefix, condition, postfix)

# plotting 
myplot_YCy(res, output_path, prefix, condition, postfix)
myplot_YcoefStd(res, output_path, prefix, condition, postfix)
myplot_cca_scatter(res, output_path, prefix, condition, postfix)
myplot_XcoefStd(res, output_path, prefix, condition, postfix)
myplot_XCx(res, output_path, prefix, condition, postfix)
