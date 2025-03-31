clc;clear
% CurrentDir = pwd;
% addpath(genpath('/oak/stanford/groups/menon/toolboxes/spm12'))

bn_atlas = '/Users/zhangyuan/Google Drive/2021_HCP_Gender_DNN/data/data_for_additional_analysis/BN_Atlas_246_2mm.nii';
bn_atlas_v=spm_vol(bn_atlas);
[bn_atlas_Y,XYZmm]=spm_read_vols(bn_atlas_v);

shier_14net = '/Users/zhangyuan/Google Drive/2021_HCP_Gender_DNN/data/data_for_additional_analysis/Shier_14Networks/*.nii';
file_list=dir(shier_14net);


map = zeros(247,15);
map(:,1) = 0:246;

cnt = zeros(247,15);
cnt(:,1) = 0:246;

for i = 1:length(file_list)
    fname = sprintf('%s/%s',file_list(i).folder,file_list(i).name);
    
    shier_v = spm_vol(fname);
    [shier_Y,XYZmm]=spm_read_vols(shier_v);
    
    bn_in_shier = bn_atlas_Y(find(shier_Y == 1));
%     file_list(i).bn_id = nonzeros(unique(bn_in_shier)); 
    
    [GC,GR] = groupcounts(bn_in_shier);
    map(GR+1,i+1) = 1;
    cnt(GR+1,i+1) = GC;
    
end

map(1,:) = [];
map(:,1) = [];
cnt(1,:) = [];
cnt(:,1) = [];

[M,I] = max(cnt,[],2);
% save('shier_to_bn.mat','file_list', 'map', 'cnt', 'M', 'I');

mat = zeros(246,2);
mat(:,1) = 1:246;
mat(:,2) = I;

net = cell(246,1);
for i = 1:246
    net{i} = file_list(I(i)).name(1:end-4);
end

fid = fopen( 'shier_to_bn.csv', 'w' );
for i = 1:246
    fprintf( fid, '%d, %d, %s\n', mat(i,1), mat(i,2), net{i});
end
fclose( fid );
