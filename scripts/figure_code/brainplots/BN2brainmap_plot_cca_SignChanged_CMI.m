% write brainbeh associations (r) to nii file and do brain plots

clc;clear;
addpath(genpath('/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/scripts/customcolormap'))
% path and files
project_path = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter';
mask_path = sprintf('%s/scripts/BN_Atlas_246_1mm.nii',project_path);
fname = sprintf('%s/results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef_cmi_n760.csv',project_path);
% fname = sprintf('%s/results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef_cmi_n760.csv',project_path);

% set up for brain plot
cfg_file = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/scripts/yz_figure_code/exp_Cfg.mat';
% preset = 'orange-white-purple'; %'jet'; %'brown-white-pool'; %'pink-white-green'; %'purple-white-green'; %'orange-white-purple'; %'red-white-blue'; %'red-yellow-green'; %'red-yellow-blue'; %'pasteljet';
preset = 'jet';

% load data
data = readtable(fname);
colnames = data.Properties.VariableNames;
ncol = length(colnames);
nroi = 218;

% load BN mask
V = spm_vol(fullfile(mask_path));
[Y,XYZmm] = spm_read_vols(V);


% for i=1:ncol
% for i=2:2:4
for i = 4

    % write to nii files
    dat = data{:, colnames{i}};
    
    Y_new = zeros(size(Y));
    
    % % without thresholding
    % for j=1:nroi
    %       % Y_new(find(Y==j)) = dat(j);
    %       Y_new(find(Y==j)) = -1 * dat(j);
    % end
    % output_file = sprintf('%s/results/cca/brainmaps/%s_ageinmodel_SignChanged_cca_cmi_n760.nii',project_path, colnames{i});

    % with thresholding
    thresh = 80;
    P = prctile(abs(dat),thresh);
    idx = find(abs(dat) > P);
    
    idx_output_f = sprintf('%s/results/cca/brainmaps/%s_thr%d_SignChanged_cca_cmi_n760.csv',project_path,colnames{i}, thresh);
    writematrix(idx(:), idx_output_f);

    for j=1:length(idx)
        % Y_new(find(Y==idx(j))) = dat(idx(j));
        Y_new(find(Y==idx(j))) = -1 * dat(idx(j));
    end
    output_file = sprintf('%s/results/cca/brainmaps/%s_thr%d_SignChanged_cca_cmi_n760.nii',project_path,colnames{i}, thresh);

    V.fname = output_file;
    V.dt = [64 0];
    V.private.dat.fname = V.fname;  
    spm_write_vol(V,Y_new); 

    % brain plot
    % update configurations
    load(cfg_file)
    EC.vol.CMt = jet(1000);
    % EC.vol.CMt = customcolormap_preset(preset);
    % orange-white-purple
    % EC.vol.CMt = customcolormap(linspace(0,1,11), {'#7f3c0a','#b35807','#e28212','#f9b967','#ffe0b2','#f7f7f5','#d7d9ee','#b3abd2','#8073a9','#562689','#2f004d'});
    EC.vol.px = ceil(max(max(dat), abs(min(dat))) * 10) / 10 + 0.1;
    EC.vol.nx = -1*EC.vol.px;
    EC.vol.pn = 0;
    EC.vol.nn = 0;
    EC.vol.pn = 1e-15;
    EC.vol.nn = -1e-15;
    save(cfg_file, 'EC')

    % plot
    nii_file = output_file;
    % output_figure = sprintf('%s/results/cca/brainmaps/%s_%s_SignChanged_cca_cmi_n760.jpg',project_path, colnames{i},preset);
    output_figure = sprintf('%s/results/cca/brainmaps/%s_%s_thr%d_SignChanged_cca_cmi_n760.jpg',project_path, colnames{i},preset,thresh);
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed_tal.nv',nii_file,cfg_file,output_figure);
end
