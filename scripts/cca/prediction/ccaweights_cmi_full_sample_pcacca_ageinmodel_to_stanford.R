setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(readxl)
library(permute)
library(reshape2)
library(ggplot2)

### function for permutation test ###
my.perm <- function(numsub, nperm, X, Y){
  
  actual_r = as.numeric(cor(X, Y))
  perm.idx = shuffleSet(numsub,nperm) 
  r = c()
  for(i in c(1:nperm)){
    idx = perm.idx[i,]
    r = c(r, cor(X[idx], Y))
  }
  p_perm = length(which(r > actual_r)) / nperm
  
  p_perm 
}
nperm = 5000

#############################################
## load Stanford brain and behavioral data ##
#############################################
# roi labels
bn_atlas_file = "data/bn_atlas.xlsx"
bn_atlas = read_excel(bn_atlas_file, sheet=1)
roi.names = bn_atlas$Description[1:218]

# behavior
beh_file = "data/subjectlists/stanford/subjectlist_baseline_wiat2_covariates_woEO.csv"
beh = read.csv(beh_file)
idx = which(beh$pid %in% c(203, 209, 7401, 7698))
beh = beh[-idx,]

# brain
brain_file = 'data/gmv_lCT_n231.csv'

brain = read.csv(brain_file)
brain = as.matrix(brain[,-c(1:3)])
sum(complete.cases(brain))
numsub = dim(brain)[1]
numroi = dim(brain)[2]


#######################################################
## load CMI Math CCA model and make prediction ##
#######################################################
# CCA results math original model
orig_path = "results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_CMI_full_math_n760/"
orig_file = "CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep=""))

# use CMI PCA rotation matrix and the mean to convert
# stanford brain data into PCs
su_pcs = sweep(brain,2,my.pca$center,'-') %*% my.pca$rotation
su_pcs = su_pcs[,c(1:my.pca$numPC)]

# apply CMI mode 2 coefficients to stanford_pcs to get the CCA mode 2 score
brain_score = su_pcs %*% res$xcoef.raw[,2] 
beh_score = beh$age * res$ycoef.raw[1,2] + beh$mathrea_std * res$ycoef.raw[2,2] +  beh$numop_std * res$ycoef.raw[3,2]

# examine correlation
cor.test(brain_score, beh_score)
p.perm = my.perm(numsub, nperm, brain_score, beh_score)
p.perm

# df = data.frame(predicted_brain_score = brain_score, predicted_beh_score = beh_score)
# write.csv(df, "figures_updated/Fig6_prediction_in_stanford/prediction_math_in_su.csv",row.names = F)


##########################################################
## load CMI Reading CCA model and make prediction ##
##########################################################
# CCA results reading original model
orig_path = "results/cca/wholebrain_rois_gmv_ageinmodel_brainnetome_CMI_full_reading_n760/"
orig_file = "CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_perm5000_pcaNoscale_ccaNoScale.RData"
load(paste(orig_path, orig_file, sep=""))

# use CMI PCA rotation matrix and the mean to convert
# Stanford brain data into PCs
su_pcs = sweep(brain,2,my.pca$center,'-') %*% my.pca$rotation
su_pcs = su_pcs[,c(1:my.pca$numPC)]

# apply CMI mode 2 coefficients to stanford_pcs to get the CCA mode 2 score
brain_score = su_pcs %*% res$xcoef.raw[,2] 
beh_score = beh$age * res$ycoef.raw[1,2] + beh$readcomp_std * res$ycoef.raw[2,2] + beh$wordread_std * res$ycoef.raw[3,2]
# examine correlation
cor.test(brain_score, beh_score)
p.perm = my.perm(numsub, nperm, brain_score, beh_score)
p.perm

# df = data.frame(predicted_brain_score = brain_score, predicted_beh_score = beh_score)
# write.csv(df, "figures_updated/Fig6_prediction_in_stanford/prediction_read_in_su.csv",row.names = F)
