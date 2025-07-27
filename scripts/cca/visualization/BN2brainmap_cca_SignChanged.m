%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Whole-Brain and Thresholded GMV Weight Maps for CCA Modes
% Author: Yuan Zhang
% Date: 2025-07-25
%
% Description:
% This script generates NIfTI brain maps and thresholded maps 
% for both math- and reading-related modes, across both CMI and Stanford cohorts. 
% It uses GMV loading weights from CCA results.
%
% For each dataset-mode combination:
%   1) Load CCA coefficients from CSV.
%   2) Map ROI weights to the Brainnetome atlas.
%   3) Save whole-brain GMV weight maps (NIfTI).
%   4) Save thresholded top-percentile GMV weight maps and their ROI indices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ----------------------------------------------------------------------
% 1. Setup Paths
% ----------------------------------------------------------------------
clc; clear;

% Project path
project_path = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub';

% Path to Brainnetome atlas mask (1mm resolution)
mask_path = sprintf('%s/data/atlas/BN_Atlas_246_1mm.nii', project_path);

% ROI count (BN atlas)
nroi = 218; 

% Output suffix (for filenames)
postfix = "ageinmodel_SignChanged_cca";

% Threshold percentile for creating top-weighted maps
thresh = 80;  % Keep weights in top 20% by absolute value

%% ----------------------------------------------------------------------
% 2. Define Datasets and Modes
% ----------------------------------------------------------------------
% Each entry in "configs" contains:
%   1) Subfolder (CMI or Stanford)
%   2) CSV file with CCA coefficients
%   3) A descriptive mode name for output files

configs = {
    % CMI Math
    'cmi', 'wholebrain_cca_cmi_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef.csv', 'cmi_math';
    
    % CMI Reading
    'cmi', 'wholebrain_cca_cmi_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef.csv', 'cmi_reading';
    
    % Stanford Math
    'stanford', 'wholebrain_cca_stanford_math/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef.csv', 'stanford_math';
    
    % Stanford Reading
    'stanford', 'wholebrain_cca_stanford_reading/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef.csv', 'stanford_reading';
};

%% ----------------------------------------------------------------------
% 3. Load Brainnetome Atlas Mask
% ----------------------------------------------------------------------
V = spm_vol(fullfile(mask_path));   % Load volume information
[Y, XYZmm] = spm_read_vols(V);      % Read voxel values (ROI indices)

%% ----------------------------------------------------------------------
% 4. Process Each Dataset-Mode
% ----------------------------------------------------------------------
for c = 1:size(configs, 1)
    dataset   = configs{c, 1};
    csv_file  = configs{c, 2};
    mode_name = configs{c, 3};
    
    fprintf('\nProcessing: %s (Dataset: %s)\n', mode_name, dataset);
    
    % Construct full path to the CCA coefficient CSV
    fname = sprintf('%s/results/cca/%s/%s', project_path, dataset, csv_file);
    
    % Load CCA coefficients table
    data = readtable(fname);
    colnames = data.Properties.VariableNames;
    
    % Find columns that correspond to brain weights for mode 2
    col_idx = find(strcmp(colnames, 'xloading_m2'));
    
    % Process each "xloading_m2" column (only 1 map, but used for loop for flexibility)
    for ii = 1:length(col_idx)
        i = col_idx(ii);
        col_name = colnames{i};
        
        fprintf('  -> Writing brain map for column: %s\n', col_name);
        
        % Extract brain weight vector for all ROIs
        dat = data{:, col_name};
        
        % ------------------------------------------------------------------
        % 4.1 Whole-Brain Map
        % ------------------------------------------------------------------
        Y_new = zeros(size(Y));
        for j = 1:nroi
            % Assign GMV weight to voxels of ROI j (inverse sign if needed)
            Y_new(Y == j) = -1 * dat(j);
        end
        
        % Output NIfTI for the whole-brain GMV weight map
        output_file = sprintf('%s/results/cca/brainmaps/%s_%s_%s.nii', ...
            project_path, col_name, mode_name, postfix);
        V.fname = output_file;
        V.dt = [64 0];  % 64 = float32
        V.private.dat.fname = V.fname;  
        spm_write_vol(V, Y_new);
        
        % ------------------------------------------------------------------
        % 4.2 Thresholded Map
        % ------------------------------------------------------------------
        % Determine threshold value (top 20% of absolute values)
        P = prctile(abs(dat), thresh);
        idx = find(abs(dat) > P);  % ROI indices above threshold
        
        % Save the indices of thresholded ROIs
        idx_output_f = sprintf('%s/results/cca/brainmaps/%s_thr%d_%s_%s.csv', ...
            project_path, col_name, thresh, mode_name, postfix);
        writematrix(idx(:), idx_output_f);
        
        % Build thresholded NIfTI map
        Y_new_thr = zeros(size(Y));
        for j = 1:length(idx)
            Y_new_thr(Y == idx(j)) = -1 * dat(idx(j));
        end
        
        % Output thresholded NIfTI
        output_file_top = sprintf('%s/results/cca/brainmaps/%s_thr%d_%s_%s.nii', ...
            project_path, col_name, thresh, mode_name, postfix);
        V.fname = output_file_top;
        V.dt = [64 0];
        V.private.dat.fname = V.fname;  
        spm_write_vol(V, Y_new_thr);
    end
end

fprintf('\nAll brain maps generated successfully!\n');
