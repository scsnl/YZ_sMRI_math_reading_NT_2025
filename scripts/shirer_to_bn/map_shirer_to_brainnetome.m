%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign Shirer Network to Brainnetome ROI
% Author: Yuan Zhang
% Date: 2025-07-25
%
% Description:
% This script maps each Brainnetome atlas ROI (BN_ROI) to a Shirer 
% 14-network parcellation based on the overlap of voxel locations. 
% For each BN ROI, the Shirer network with the maximum voxel count 
% is assigned.
%
% Outputs:
%   - shirer_to_bn.mat: Contains mapping matrices and summary info.
%   - shirer_to_bn.csv: CSV file listing each BN ROI and its assigned network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -----------------------------------------------------------------
% 1. Setup
% ------------------------------------------------------------------

clc; clear;

% Path to the project folder
project_path = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub/';

% ------------------------------------------------------------------
% 2. Load Brainnetome Atlas
% ------------------------------------------------------------------

% Load BN atlas (NIfTI file)
bn_atlas = sprintf('%s/data/atlas/BN_Atlas_246_2mm.nii', project_path);
bn_atlas_v = spm_vol(bn_atlas);             % Get volume information
[bn_atlas_Y, XYZmm] = spm_read_vols(bn_atlas_v);  % Read voxel data

% ------------------------------------------------------------------
% 3. Load Shirer 14-Network Atlas Files
% ------------------------------------------------------------------

% Get all Shirer 14-network atlas NIfTI files
shirer_14net = sprintf('%s/data/atlas/Shirer_14Networks/*.nii', project_path);
file_list = dir(shirer_14net);  % Structure array of Shirer atlas files

% ------------------------------------------------------------------
% 4. Initialize Mapping Matrices
% ------------------------------------------------------------------
% map: Binary matrix indicating if BN ROI (row) belongs to a Shirer network (column)
% cnt: Count matrix showing number of voxels of BN ROI that overlap with Shirer network

map = zeros(247, 15);   % 246 BN ROIs + 1 header row
map(:,1) = 0:246;       % First column: ROI indices (0 to 246)

cnt = zeros(247, 15);
cnt(:,1) = 0:246;

% ------------------------------------------------------------------
% 5. Map BN ROIs to Shirer Networks
% ------------------------------------------------------------------

for i = 1:length(file_list)
    % Load the Shirer network mask
    fname = sprintf('%s/%s', file_list(i).folder, file_list(i).name);
    shirer_v = spm_vol(fname);
    [shirer_Y, XYZmm] = spm_read_vols(shirer_v);
    
    % Find BN ROIs that overlap with the current Shirer network
    bn_in_shirer = bn_atlas_Y(find(shirer_Y == 1));
    file_list(i).bn_id = nonzeros(unique(bn_in_shirer));  % Store ROI IDs
    
    % Count the number of voxels for each BN ROI in this Shirer network
    [GC, GR] = groupcounts(bn_in_shirer);
    
    % Mark overlap (map) and store voxel count (cnt)
    map(GR + 1, i + 1) = 1;
    cnt(GR + 1, i + 1) = GC;
end

% Remove the "0" ROI row and column (background)
map(1,:) = [];
map(:,1) = [];
cnt(1,:) = [];
cnt(:,1) = [];

% ------------------------------------------------------------------
% 6. Determine Shirer Network Assignment
% ------------------------------------------------------------------

% For each BN ROI, find the Shirer network with maximum voxel count
[M, I] = max(cnt, [], 2);

% Save the raw mapping data
save(sprintf('%s/data/atlas/shirer_to_bn.mat', project_path), 'file_list', 'map', 'cnt', 'M', 'I');

% Create matrix (BN ROI index vs assigned Shirer network index)
mat = zeros(246, 2);
mat(:,1) = 1:246;
mat(:,2) = I;

% Create a cell array of assigned Shirer network names
net = cell(246, 1);
for i = 1:246
    net{i} = file_list(I(i)).name(1:end-4);  % Remove ".nii" extension
end

% ------------------------------------------------------------------
% 7. Save Results to CSV
% ------------------------------------------------------------------

fid = fopen(sprintf('%s/data/atlas/shirer_to_bn.csv', project_path), 'w');
for i = 1:246
    fprintf(fid, '%d, %d, %s\n', mat(i,1), mat(i,2), net{i});
end
fclose(fid);
