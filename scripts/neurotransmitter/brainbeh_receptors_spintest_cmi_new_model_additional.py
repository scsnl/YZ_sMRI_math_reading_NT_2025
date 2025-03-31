
from utility import *
import pandas as pd
# from matplotlib.colors import ListedColormap
from scipy.stats import zscore

"""
Path and files
"""

# path = '/Users/zhangyuan/Google Drive/2020_LongitudinalStructure_Reading_Math/hansen_receptors-main/'
project_dir = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/'
brainnetome = '/Users/zhangyuan/Google Drive/2020_LongitudinalStructure_Reading_Math/BN_Atlas_246_2mm.nii.gz'
scale = 'bn246'
dist_fname = project_dir + 'scripts/spin/BNA_distmat_218.txt'
n_perm = 5000  # number of permutations for null models
# spatnull = 'burt2020'
spatnull = 'moran'

"""
load brain-behavioral association maps
"""

# CCA brain-behavior map
fname = project_dir + '/results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_coef_cmi_n760.csv'
# fname = project_dir + '/results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_coef_cmi_n760.csv'
output_file = project_dir + '/results/cca/CCA_PCA_roi_gmv_brainnetome_mathstd_ageinmodel_individual_neurotransmitter_regression_results_cmi_n760.csv'
# output_file = project_dir + '/results/cca/CCA_PCA_roi_gmv_brainnetome_readstd_ageinmodel_individual_neurotransmitter_regression_results_cmi_n760.csv'
# idx_list = [1,3]
idx_list = [3]

bmaps = pd.read_csv(fname)
bmap_names = list(bmaps.columns.values)
bmaps = bmaps.to_numpy()

"""
load receptor data
"""
receptor_data_complete = np.genfromtxt(project_dir +'scripts/neurotransmitter/receptor_data_'+scale+'.csv', delimiter=',')[:218,]
receptor_names = np.load(project_dir +'scripts/neurotransmitter/receptor_names_pet_yz.npy')
# Store different sets of receptors of interest in a dictionary
# receptor_sets = {
#     'dopamine': ['D1', 'D2', 'DAT'],
#     'glutamate': ['mGluR5', 'NMDA'],
#     'GABA': ['GABAa'],
#     'serotonin': ['5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT'],
#     'acetylcholine': ['A4B2','M1','VAChT'],
#     'norepinephrine': ['NET'],
#     'cannabinoid': ['CB1'],
#     'histamine': ['H3'],
#     'opioid': ['MOR'],
#     'DGG': ['D1', 'D2', 'DAT', 'GABAa', 'mGluR5', 'NMDA'],
#     'all_13': ['5HT1a','5HT1b','5HT2a','5HT4','5HT6','5HTT','A4B2','M1','VAChT','NET','CB1','H3','MOR'],
#     'all_19': ['5HT1a','5HT1b','5HT2a','5HT4','5HT6','5HTT','A4B2','CB1','D1', 'D2', 'DAT','GABAa',
#                'H3','M1','mGluR5','MOR','NET','NMDA','VAChT']
# }

receptor_sets = {
    'D1': ['D1'],
    'D2': ['D2'],
    'DAT': ['DAT'],
    'mGluR5': ['mGluR5'],
    'NMDA': ['NMDA'],
    'GABAa': ['GABAa'],
    '5HT1a': ['5HT1a'],
    '5HT1b': ['5HT1b'],
    '5HT2a': ['5HT2a'],
    '5HT4': ['5HT4'],
    '5HT6': ['5HT6'],
    '5HTT': ['5HTT'],
    'A4B2': ['A4B2'],
    'M1': ['M1'],
    'VAChT': ['VAChT'],
    'NET': ['NET'],
    'CB1': ['CB1'],
    'H3': ['H3'],
    'MOR': ['MOR']
}

"""
multiple linear regression with spin test
"""
# Initialize lists to store results
results = {
    'receptor_set': [],
    'brain_map': [],
    'adj_r2': [],
    'moran_p': []
}

# Iterate through the dictionary and conduct analysis
# count = 0
for set_name, receptors in receptor_sets.items():
    print(f"Analyzing set: {set_name}")
    print(f"Receptors: {receptors}")

    # get the subset of receptor data
    idx = np.where(np.in1d(receptor_names, receptors))[0]
    receptor_data = receptor_data_complete[:, idx]
    print(receptor_data.shape)
    if receptor_data.ndim == 1:
        receptor_data = receptor_data.reshape(-1, 1)

    for i in range(len(idx_list)):
        print(idx_list[i])
        print(bmap_names[idx_list[i]])

        # get adjust r2 of the model
        r2, coef = get_reg_r_sq(zscore(receptor_data), zscore(-1*bmaps[:, idx_list[i]]))
        print('adjusted r square: {}'.format(r2))
        print('regression coefficients: {}'.format(coef))

        # get model pval (spin test)
        surrogates = make_surrogates(bmaps[:,idx_list[i]], dist_fname, spatnull, n_perm) # (N, `N_PERM`) np.ndarray
        print(surrogates.shape)
        pval = get_reg_r_pval(zscore(receptor_data), zscore(bmaps[:, idx_list[i]]), zscore(surrogates))
        print('p val (spin test): {}'.format(pval))

        # Store the results in the lists
        results['receptor_set'].append(set_name)
        results['brain_map'].append(bmap_names[idx_list[i]])
        results['adj_r2'].append(r2)
        results['moran_p'].append(pval)

    # count += 1

# Convert the results dictionary to a DataFrame
results_df = pd.DataFrame(results)
# Save the DataFrame to a CSV file
results_df.to_csv(output_file, index=False)
print("Results saved")