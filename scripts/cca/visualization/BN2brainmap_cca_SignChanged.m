%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Whole-Brain GMV Weight Maps for CCA Modes
% Author: Yuan Zhang
% Date: 2025-07-25
%
% Description:
% This script generates whole-brain NIfTI brain maps for both math- and
% reading-related modes, across both CMI and Stanford cohorts.
% It uses GMV loading weights from CCA results.
%
% For each dataset-mode combination:
%   1) Load CCA coefficients from CSV.
%   2) Map ROI weights to the Brainnetome atlas.
%   3) Save whole-brain GMV weight maps as float32 compressed NIfTI files.
%   4) Save ROI-level source data used to generate the NIfTI maps.
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

% Output directory
brainmap_dir = sprintf('%s/results/cca/brainmaps', project_path);
if ~exist(brainmap_dir, 'dir')
    mkdir(brainmap_dir);
end

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
    
    % Process each "xloading_m2" column
    for ii = 1:length(col_idx)
        i = col_idx(ii);
        col_name = colnames{i};
        
        fprintf('  -> Writing brain map for column: %s\n', col_name);
        
        % Extract brain weight vector for all ROIs
        dat = data{:, col_name};
        
        % Sign-oriented ROI values actually used for NIfTI visualization
        dat_sign_oriented = -1 * dat;
        
        % ------------------------------------------------------------------
        % 4.1 Save ROI-level source data used for NIfTI map
        % ------------------------------------------------------------------
        roi_id = (1:nroi)';
        
        source_table = table( ...
            roi_id, ...
            dat(:), ...
            dat_sign_oriented(:), ...
            'VariableNames', { ...
                'ROI_ID', ...
                'Original_xloading_m2', ...
                'Sign_oriented_value_used_for_NIfTI' ...
            } ...
        );
        
        source_csv_file = sprintf('%s/%s_%s_%s_source_data_roi_values.csv', ...
            brainmap_dir, col_name, mode_name, postfix);
        
        writetable(source_table, source_csv_file);
        
        % ------------------------------------------------------------------
        % 4.2 Whole-Brain Map
        % ------------------------------------------------------------------
        Y_new = zeros(size(Y), 'single');  % float32 array
        
        for j = 1:nroi
            % Assign sign-oriented GMV weight to voxels of ROI j
            Y_new(Y == j) = single(dat_sign_oriented(j));
        end
        
        % Output NIfTI for the whole-brain GMV weight map
        output_file = sprintf('%s/%s_%s_%s.nii', ...
            brainmap_dir, col_name, mode_name, postfix);
        
        Vout = V;
        Vout.fname = output_file;
        Vout.dt = [16 0];  % 16 = float32 in SPM
        Vout.private.dat.fname = Vout.fname;
        
        spm_write_vol(Vout, Y_new);
        
        % Compress to .nii.gz and remove uncompressed .nii
        gzip(output_file);
        delete(output_file);
    end
end

fprintf('\nAll brain maps generated successfully!\n');