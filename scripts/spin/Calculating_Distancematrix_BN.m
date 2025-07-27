
% -------------------------------------------------------------------------
% Brainnetome Atlas Region Coordinates and Distance Matrix
% Author: Yuan Zhang
% Date: 2025-07-25
%
% Description:
% This script:
%   1. Loads Brainnetome (BN) atlas subregion information (IDs and MNI coordinates).
%   2. Extracts and combines left and right hemisphere region indices and coordinates.
%   3. Creates a matrix of all region IDs (1–218) with corresponding MNI coordinates.
%   4. Saves the region coordinates to a text file.
%   5. Computes the pairwise Euclidean distance matrix between all regions.
%   6. Saves the distance matrix to a text file.
%
% Outputs:
%   - BNA_coordinates_218.txt: 218 x 3 matrix of MNI coordinates.
%   - BNA_distmat_218.txt: 218 x 218 matrix of Euclidean distances between regions.
% -------------------------------------------------------------------------

clear,clc

% -------------------------------------------------------------------------
% 1. Project setup
% -------------------------------------------------------------------------

project_dir = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub';

% Load Brainnetome atlas subregions table
T = readtable(sprintf('%s/data/atlas/BNA_subregions.xlsx', project_dir));

% -------------------------------------------------------------------------
% 2. Define region indices
% -------------------------------------------------------------------------

num_regions = 218; % Total number of BN regions considered
ee = num_regions/2; % Number of regions per hemisphere (109 left + 109 right)

% Extract region indices for left and right hemispheres
region_ind_L = T.LabelID_L(1:ee);
region_ind_R = T.LabelID_R(1:ee);

% -------------------------------------------------------------------------
% 3. Parse MNI coordinates
% -------------------------------------------------------------------------

% Convert left hemisphere MNI coordinates from text to numeric array
region_MNI_L = cell2mat(cellfun(@str2num,T.lh_MNI_X_Y_Z_(1:ee),'uniform',0));

% Convert right hemisphere MNI coordinates from text to numeric array
region_MNI_R = cell2mat(cellfun(@str2num,T.rh_MNI_X_Y_Z_(1:ee),'uniform',0));

% -------------------------------------------------------------------------
% 4. Combine hemispheres
% -------------------------------------------------------------------------

region_ind_all=[region_ind_L;region_ind_R]; % Combine region indices
region_MNI_all=[region_MNI_L;region_MNI_R]; % Combine coordinates

% -------------------------------------------------------------------------
% 5. Create final region index and coordinate matrix
% -------------------------------------------------------------------------
% region_final: 
%   Column 1 = Region ID (1 to 218)
%   Columns 2-4 = MNI coordinates (X, Y, Z)

for i=1:num_regions
   ind_temp=find(region_ind_all==i); % Find corresponding index
   region_final(i,1)=i; % Assign region ID
   region_final(i,2:4)=region_MNI_all(ind_temp,:); % Assign MNI coordinates
end

% Extract only the coordinates (X, Y, Z)
region_coord = region_final(:,2:4);

% Save coordinates to text file
save(sprintf('%s/scripts/spin/BNA_coordinates_218.txt', project_dir),'region_coord','-ascii');

% -------------------------------------------------------------------------
% 6. Compute pairwise distance matrix
% -------------------------------------------------------------------------
% D(j, k) = Euclidean distance between region j and region k

for j=1:length(region_final)
    for k=1:length(region_final)
        D(j,k) = pdist2(region_final(j,2:4),region_final(k,2:4));
    end
end

% Save distance matrix to text file
save(sprintf('%s/scripts/spin/BNA_distmat_218.txt', project_dir),'D','-ascii');
