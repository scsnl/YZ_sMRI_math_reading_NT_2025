# -*- coding: utf-8 -*-
"""
Concatenate parcellated PET images into a region x receptor matrix of densities.
Author: Yuan Zhang
Date: 2025-07-25

Steps:
1. Load parcellated PET receptor data (already reduced to BN246 regions).
2. Combine these data into a matrix of shape (n_regions x n_receptors).
3. Average datasets for the same receptor type, weighting by sample size (e.g., 5HT1B).
4. Save receptor names and final receptor data for further analyses.
"""

import numpy as np
from scipy.stats import zscore

# -------------------------------------------------------------------------
# 1. Setup
# -------------------------------------------------------------------------

# Name of the parcellation scale (Brainnetome 246 regions)
scale = 'bn246'

# Project base directory
path = '/Users/zhangyuan/Google Drive/2020_LongitudinalStructure_Reading_Math/GitHub/'

# Number of nodes (BN246 regions)
nnodes = 246

# -------------------------------------------------------------------------
# 2. List of PET Receptor CSVs
# -------------------------------------------------------------------------

# These CSV files contain the parcellated PET values (one column per region).
receptors_csv = [path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HT1a_way_hc36_savli.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HT1b_p943_hc22_savli.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HT1b_p943_hc65_gallezot.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HT2a_cimbi_hc29_beliveau.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HT4_sb20_hc59_beliveau.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HT6_gsk_hc30_radhakrishnan.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/5HTT_dasb_hc100_beliveau.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/A4B2_flubatine_hc30_hillmer.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/CB1_omar_hc77_normandin.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/D1_SCH23390_hc13_kaller.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/D2_flb457_hc37_smith.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/D2_flb457_hc55_sandiego.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/DAT_fpcit_hc174_dukart_spect.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/GABAa-bz_flumazenil_hc16_norgaard.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/H3_cban_hc8_gallezot.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/M1_lsn_hc24_naganawa.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/mGluR5_abp_hc22_rosaneto.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/mGluR5_abp_hc28_dubois.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/mGluR5_abp_hc73_smart.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/MU_carfentanil_hc204_kantonen.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/NAT_MRB_hc77_ding.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/NMDA_ge179_hc29_galovic.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/VAChT_feobv_hc4_tuominen.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/VAChT_feobv_hc5_bedard_sum.csv',
                 path+'data/neurotransmitter/PET_parcellated/'+scale+'/VAChT_feobv_hc18_aghourian_sum.csv']

# -------------------------------------------------------------------------
# 3. Combine All Receptor Data
# -------------------------------------------------------------------------
# Initialize matrix (n_regions x n_receptor_files)
r = np.zeros([nnodes, len(receptors_csv)])
print("Initial receptor matrix shape:", r.shape)

# Fill the matrix with values from each receptor CSV
for i in range(len(receptors_csv)):
    r[:, i] = np.genfromtxt(receptors_csv[i], delimiter=',')

# Define unique receptor names (grouping by receptor type)
receptor_names = np.array(["5HT1a", "5HT1b", "5HT2a", "5HT4", "5HT6", "5HTT", "A4B2",
                           "CB1", "D1", "D2", "DAT", "GABAa", "H3", "M1", "mGluR5",
                           "MOR", "NET", "NMDA", "VAChT"])

# Save receptor names for later use
np.save(path+'data/neurotransmitter/receptor_names_pet.npy', receptor_names)

# -------------------------------------------------------------------------
# 4. Construct Final Region x Receptor Matrix
# -------------------------------------------------------------------------

# Initialize final matrix with shape (n_regions x n_receptors)
receptor_data = np.zeros([nnodes, len(receptor_names)])

# Map columns from raw data (r) to receptor_data based on receptor type
receptor_data[:, 0] = r[:, 0]        # 5HT1a
receptor_data[:, 2:9] = r[:, 3:10]   # 5HT2a through CB1
receptor_data[:, 10:14] = r[:, 12:16] # DAT through M1
receptor_data[:, 15:18] = r[:, 19:22] # MOR through NMDA

# -------------------------------------------------------------------------
# 5. Weighted Averages for Receptors with Multiple Datasets
# -------------------------------------------------------------------------
# Weighted average of 5HT1B (p943) combining 2 datasets (22 & 65 subjects)
receptor_data[:, 1] = (zscore(r[:, 1])*22 + zscore(r[:, 2])*65) / (22+65)

# Weighted average of D2 (flb457) combining 2 datasets (37 & 55 subjects)
receptor_data[:, 9] = (zscore(r[:, 10])*37 + zscore(r[:, 11])*55) / (37+55)

# Weighted average of mGluR5 (ABP688) combining 3 datasets (22, 28, 73 subjects)
receptor_data[:, 14] = (zscore(r[:, 16])*22 + zscore(r[:, 17])*28 + zscore(r[:, 18])*73) / (22+28+73)

# Weighted average of VAChT (FEOBV) combining 3 datasets (4, 5, 18 subjects)
receptor_data[:, 18] = (zscore(r[:, 22])*4 + zscore(r[:, 23]) + zscore(r[:, 24])) / (4+5+18)

# -------------------------------------------------------------------------
# 6. Save Final Matrix
# -------------------------------------------------------------------------

np.savetxt(path+'data/neurotransmitter/receptor_data_'+scale+'.csv', receptor_data, delimiter=',')
print("Receptor data matrix saved successfully.")