setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter")
rm(list=ls())

library(stringr)
library(readxl)

# path
data_path = 'data/lCT/CMI/summary_stats_output_cmihbn_expand'
bn_atlas_file = "data/bn_atlas.xlsx"
beh_file = 'data/subjectlists/cmi/2024_04_22_CMI_supercleanTD_completeNP_7-13y_n858_withNPinNeed_updatename.csv'

# load data
# roi labels
bn_atlas = read_excel(bn_atlas_file, sheet=1)
roi.names = bn_atlas$Description[1:218]
#
beh = read.csv(beh_file)
nsub = dim(beh)[1]
count = 0
missing_pids = c()
existing_pids = c()
good_pids = c()
bad_pids = c()

for(i in c(1:nsub)){
  fname = paste(data_path,'/',beh$oak_id[i], 
                '_visit0_session0_246_BN_roi_1mm_ct_stats_wdelim.csv',sep="")
  if(file.exists(fname)){
    indiv = read.csv(fname, header=TRUE)
    col_of_interest = c("Label", "Mean", "SurfaceAreaInMillimetersSquared",
                        "VolumeInVoxels", "dimensions_x", "dimensions_y", "dimensions_z")
    idx = which(colnames(indiv) %in% col_of_interest)
    indiv = indiv[,idx]
    indiv = indiv[!duplicated(indiv), ]
    
    existing_pids = c(existing_pids, beh$oak_id[i])
    
    indiv_ct = indiv$Mean
    indiv_sa = indiv$SurfaceAreaInMillimetersSquared
    indiv_gmv = indiv$VolumeInVoxels*(indiv$dimensions_x*indiv$dimensions_y*indiv$dimensions_z)
    
    nroi_ct = dim(indiv)[1]
    if(nroi_ct < 246){
       print(paste('too many ROIs are missing for ', beh$oak_id[i], ' visit0 session0', sep=""))
       bad_pids = c(bad_pids, beh$oak_id[i])
       next
    }
    
    if(sum(indiv_ct[1:218] < 0.1) > 10){
       print(paste('too many ROIs with near-0 ct for ', beh$oak_id[i], ' visit0 session0', sep=""))
       bad_pids = c(bad_pids, beh$oak_id[i])
    }else{
       good_pids = c(good_pids, beh$oak_id[i])
       if(count == 0){
          ct_all = indiv_ct
          sa_all = indiv_sa
          gmv_all = indiv_gmv
       }else{
         tryCatch({
           ct_all = rbind(ct_all, indiv_ct)
         }, warning = function(w) {
           cat(sprintf("Warning for CT at iteration %d for %s: %s\n", i, beh$oak_id[i], w$message))
         })
          
         tryCatch({
           sa_all = rbind(sa_all, indiv_sa)
         }, warning = function(w) {
           cat(sprintf("Warning for SA at iteration %d for %s: %s\n", i, beh$oak_id[i], w$message))
         })
         
         tryCatch({
           gmv_all = rbind(gmv_all, indiv_gmv)
         }, warning = function(w) {
           cat(sprintf("Warning for GMV at iteration %d for %s: %s\n", i, beh$oak_id[i], w$message))
         })

       }
       count = count + 1
    }
 
  }else{
    print(paste('data not exist for ', beh$oak_id[i], ' visit0 session0', sep=""))
    missing_pids = c(missing_pids, beh$oak_id[i])
  }
}

write.table(missing_pids, 'missing_pids_83.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(existing_pids, 'existing_pids_775.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(existing_pids, 'good_pids_760.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)

ct_all = as.data.frame(ct_all)
colnames(ct_all) = bn_atlas$`Region Label`
rownames(ct_all) = c()
ct_all = ct_all[,-c(219:246)]
ct_all$pid = good_pids
ct_all$visit = 0
ct_all$session = 0
ct_all = ct_all[,c(219:221,1:218)]

sa_all = as.data.frame(sa_all)
colnames(sa_all) = bn_atlas$`Region Label`
rownames(sa_all) = c()
sa_all = sa_all[,-c(219:246)]
sa_all$pid = good_pids
sa_all$visit = 0
sa_all$session = 0
sa_all = sa_all[,c(219:221,1:218)]

gmv_all = as.data.frame(gmv_all)
colnames(gmv_all) = bn_atlas$`Region Label`
rownames(gmv_all) = c()
gmv_all = gmv_all[,-c(219:246)]
gmv_all$pid = good_pids
gmv_all$visit = 0
gmv_all$session = 0
gmv_all = gmv_all[,c(219:221,1:218)]


outputfile = 'data/ct_lCT_cmi_n760.csv'
write.table(ct_all, outputfile, sep=",", row.names=FALSE, quote = FALSE)

outputfile = 'data/sa_lCT_cmi_n760.csv'
write.table(sa_all, outputfile, sep=",", row.names=FALSE, quote = FALSE)

outputfile = 'data/gmv_lCT_cmi_n760.csv'
write.table(gmv_all, outputfile, sep=",", row.names=FALSE, quote = FALSE)

