setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(stringr)
library(readxl)

# path
data_path = 'data/lCT/Stanford/wiat2/summary_stats_output_wiat2'
bn_atlas_file = "data/bn_atlas.xlsx"
beh_file = "data/subjectlists/stanford/subjectlist_baseline_wiat2_covariates_woEO.csv"

# load data
# roi labels
bn_atlas = read_excel(bn_atlas_file, sheet=1)
roi.names = bn_atlas$Description[1:218]
#
beh = read.csv(beh_file)
nsub = dim(beh)[1]

for(i in c(1:nsub)){
  
  fname = paste(data_path,'/',str_pad(beh$pid[i],4,pad="0"), 
                '_visit', beh$visit[i], 
                '_session', beh$session[i],
                '_246_BN_roi_1mm_ct_stats_wdelim.csv',sep="")
  if(file.exists(fname)){
    indiv = read.csv(fname, header=TRUE)
    if(i == 1){
      ct_all = c(beh$pid[i], beh$visit[i], beh$session[i],indiv$Mean) 
      sa_all = c(beh$pid[i], beh$visit[i], beh$session[i],indiv$SurfaceAreaInMillimetersSquared) 
      gmv_all = c(beh$pid[i], beh$visit[i], beh$session[i],indiv$VolumeInVoxels*(indiv$dimensions_x*indiv$dimensions_y*indiv$dimensions_z)) 
    }else{
      ct_all = rbind(ct_all, c(beh$pid[i], beh$visit[i], beh$session[i],indiv$Mean))
      sa_all = rbind(sa_all, c(beh$pid[i], beh$visit[i], beh$session[i],indiv$SurfaceAreaInMillimetersSquared))
      gmv_all = rbind(gmv_all, c(beh$pid[i], beh$visit[i], beh$session[i],indiv$VolumeInVoxels*(indiv$dimensions_x*indiv$dimensions_y*indiv$dimensions_z)))
    }
  }else{
    print(paste('data not exist for pid', beh$pid[i], ' visit', beh$visit[i], ' session', beh$session[i], sep=""))
  }
}

ct_all = as.data.frame(ct_all)
colnames(ct_all) = c("pid","visit","session",bn_atlas$`Region Label`)
rownames(ct_all) = c()
ct_all = ct_all[,-c(222:249)]

sa_all = as.data.frame(sa_all)
colnames(sa_all) = c("pid","visit","session",bn_atlas$`Region Label`)
rownames(sa_all) = c()
sa_all = sa_all[,-c(222:249)]

gmv_all = as.data.frame(gmv_all)
colnames(gmv_all) = c("pid","visit","session",bn_atlas$`Region Label`)
rownames(gmv_all) = c()
gmv_all = gmv_all[,-c(222:249)]

# Among the original subjectlist of 235 subjects, 
# 209_visit1_session1 LCT CT is empty, SA is missing
# 203_visit1_session1, 7401_visit2_session1, 7698_visit1_session1 
# did not survive ANTs pipeline

idx = which(gmv_all$pid %in% c(203, 209, 7401, 7698))
ct_all = ct_all[-idx,]
sa_all = sa_all[-idx,]
gmv_all = gmv_all[-idx,]

outputfile = 'data/ct_lCT_stanford_n231.csv'
write.table(ct_all, outputfile, sep=",", row.names=FALSE, quote = FALSE)

outputfile = 'data/sa_lCT_stanford_n231.csv'
write.table(sa_all, outputfile, sep=",", row.names=FALSE, quote = FALSE)

outputfile = 'data/gmv_lCT_stanford_n231.csv'
write.table(gmv_all, outputfile, sep=",", row.names=FALSE, quote = FALSE)


