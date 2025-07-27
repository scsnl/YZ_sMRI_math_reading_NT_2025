from utility import *
import pandas as pd
from scipy.stats import zscore
import numpy as np
import os

"""
Author: Yuan Zhang
Date: 2025-07-25

This script performs receptor-based regression analysis for:
  1. CMI math
  2. CMI reading
  3. Stanford math
  4. Stanford reading

It uses:
  - GMV weight maps from CCA.
  - PET-based receptor density maps.
  - Spin-test (spatial permutation) for significance testing.

The results for each dataset are saved as separate CSV files.
"""

# -------------------------------------------------------------------------
# Paths and general settings
# -------------------------------------------------------------------------
project_dir = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub/'
brainnetome = project_dir + 'data/atlas/BN_Atlas_246_2mm.nii.gz'
scale = 'bn246'
dist_fname = project_dir + 'scripts/spin/BNA_distmat_218.txt'
n_perm = 5000  # number of permutations for null models
spatnull = 'moran'  # null model type

# Index list: 0 = full sample, 1 = 7–8.9 years, 2 = 9–10.9 years, 3 = 11–13.9 years
idx_list = [0, 1, 2, 3]

# -------------------------------------------------------------------------
# Load receptor data
# -------------------------------------------------------------------------
receptor_data_complete = np.genfromtxt(
    project_dir + f'data/neurotransmitter/receptor_data_{scale}.csv',
    delimiter=','
)[:218,]
receptor_names = np.load(project_dir + f'data/neurotransmitter/receptor_names_pet.npy')

# Define receptor sets
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

# -------------------------------------------------------------------------
# Function to run analysis for a single dataset
# -------------------------------------------------------------------------
def run_analysis(brain_map_file, output_file, dataset_name):
    """
    Perform receptor regression analysis with spin test for one dataset.

    Parameters
    ----------
    brain_map_file : str
        Path to CSV file containing GMV weight maps for full sample and age bins.
    output_file : str
        Path where the results CSV will be saved.
    dataset_name : str
        Name of the dataset (e.g., 'CMI-math') for logging purposes.
    """
    print(f"\n=== Running analysis for {dataset_name} ===")

    # Load brain-behavior association maps
    bmaps = pd.read_csv(brain_map_file)
    bmap_names = list(bmaps.columns.values)
    bmaps = bmaps.to_numpy()

    # Prepare results dictionary
    results = {
        'receptor_set': [],
        'brain_map': [],
        'adj_r2': [],
        'moran_p': [],
        'regression_coef': []
    }

    # Iterate through each receptor set
    for set_name, receptors in receptor_sets.items():
        print(f"Analyzing receptor set: {set_name} -> {receptors}")

        # Extract receptor data subset
        idx = np.where(np.in1d(receptor_names, receptors))[0]
        receptor_data = receptor_data_complete[:, idx]

        if receptor_data.ndim == 1:
            receptor_data = receptor_data.reshape(-1, 1)

        # Loop over each age bin / brain map
        for i in idx_list:
            print(f"  Brain map: {bmap_names[i]}")

            # Regression (adjusted R2 and coefficients)
            r2, coef = get_reg_r_sq(zscore(receptor_data), zscore(bmaps[:, i]))
            print(f'    Adjusted R²: {r2:.4f}')
            print(f'    Coefficients: {coef}')

            # Spin-test (permutation-based p-value)
            surrogates = make_surrogates(bmaps[:, i], dist_fname, spatnull, n_perm)
            pval = get_reg_r_pval(zscore(receptor_data), zscore(bmaps[:, i]), zscore(surrogates))
            print(f'    Spin-test p-value: {pval:.4f}')

            # Store results
            results['receptor_set'].append(set_name)
            results['brain_map'].append(bmap_names[i])
            results['adj_r2'].append(r2)
            results['moran_p'].append(pval)
            results['regression_coef'].append(coef.item() if coef.size == 1 else coef.tolist())

    # Convert to DataFrame and save
    results_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# -------------------------------------------------------------------------
#  Run analysis for all datasets (CMI & Stanford, math & reading)
# -------------------------------------------------------------------------
datasets = [
    {
        "brain_map_file": project_dir + 'results/age_analysis/cmi/gmv_weights_agebin_math.csv',
        "output_file": project_dir + 'results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_math.csv',
        "name": "CMI-math"
    },
    {
        "brain_map_file": project_dir + 'results/age_analysis/cmi/gmv_weights_agebin_read.csv',
        "output_file": project_dir + 'results/neurotransmitter/cmi/individual_neurotransmitter_regression_results_agebin_read.csv',
        "name": "CMI-reading"
    },
    {
        "brain_map_file": project_dir + 'results/age_analysis/stanford/gmv_weights_agebin_math.csv',
        "output_file": project_dir + 'results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_math.csv',
        "name": "Stanford-math"
    },
    {
        "brain_map_file": project_dir + 'results/age_analysis/stanford/gmv_weights_agebin_read.csv',
        "output_file": project_dir + 'results/neurotransmitter/stanford/individual_neurotransmitter_regression_results_agebin_read.csv',
        "name": "Stanford-reading"
    }
]

# Execute analysis for each dataset
for ds in datasets:
    run_analysis(ds["brain_map_file"], ds["output_file"], ds["name"])
