clear,clc

T = readtable('/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/scripts/spin/BNA_subregions.xlsx')
num_regions = 218;
ee = num_regions/2;
region_ind_L = T.LabelID_L(1:ee);
region_ind_R = T.LabelID_R(1:ee);

region_MNI_L = cell2mat(cellfun(@str2num,T.lh_MNI_X_Y_Z_(1:ee),'uniform',0));
region_MNI_R = cell2mat(cellfun(@str2num,T.rh_MNI_X_Y_Z_(1:ee),'uniform',0));

region_ind_all=[region_ind_L;region_ind_R];
region_MNI_all=[region_MNI_L;region_MNI_R];
for i=1:num_regions
   ind_temp=find(region_ind_all==i)
   region_final(i,1)=i;
   region_final(i,2:4)=region_MNI_all(ind_temp,:) 
end

region_coord = region_final(:,2:4);
save('/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/scripts/spin/BNA_coordinates_218.txt','region_coord','-ascii')

for j=1:length(region_final)
    for k=1:length(region_final)
        D(j,k) = pdist2(region_final(j,2:4),region_final(k,2:4));
    end
end
save('/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/scripts/spin/BNA_distmat_218.txt','D','-ascii')
