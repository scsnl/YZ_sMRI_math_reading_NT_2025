from utility import *
import pandas as pd
from scipy.stats import zscore
import numpy as np
import os

"""
Author: Yuan Zhang
Date: 2026-04-13

This script performs receptor-based regression analysis using
the Mode 2 GMV weight map extracted from the joint CCA model
(math + reading combined model).

Specifically:
  - load *_coef.csv from the joint CCA result folder
  - extract xloading_m2
  - multiply by -1 to flip the sign
  - run receptor regression and spin-test

The results are saved separately for:
  1. CMI joint CCA shared mode (Mode 2)
  2. Stanford joint CCA shared mode (Mode 2)
"""

# -------------------------------------------------------------------------
# Paths and general settings
# -------------------------------------------------------------------------
project_dir = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub/'
scale = 'bn246'
dist_fname = project_dir + 'scripts/spin/BNA_distmat_218.txt'
n_perm = 5000
spatnull = 'moran'

# -------------------------------------------------------------------------
# Load receptor data
# -------------------------------------------------------------------------
receptor_data_complete = np.genfromtxt(
    project_dir + f'data/neurotransmitter/receptor_data_{scale}.csv',
    delimiter=','
)[:218, :]
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
# Function to run analysis for a single Mode 2 brain map
# -------------------------------------------------------------------------
def run_analysis_from_coef_csv(coef_file, output_file, dataset_name,
                               weight_col='xloading_m2', flip_sign=True):
    """
    Perform receptor regression analysis with spin test for one dataset.

    Parameters
    ----------
    coef_file : str
        Path to CCA coefficient CSV file containing xloading_m2.
    output_file : str
        Path where the results CSV will be saved.
    dataset_name : str
        Name of the dataset for logging.
    weight_col : str
        Column name of the desired brain weight map. Default = 'xloading_m2'.
    flip_sign : bool
        Whether to multiply the extracted map by -1. Default = True.
    """
    print(f"\n=== Running analysis for {dataset_name} ===")

    # ------------------------------------------------------------------
    # Load coefficient file and extract brain map
    # ------------------------------------------------------------------
    coef_df = pd.read_csv(coef_file)

    if weight_col not in coef_df.columns:
        raise ValueError(f"{weight_col} not found in {coef_file}")

    bmap = coef_df[weight_col].to_numpy(dtype=float)

    if flip_sign:
        bmap = -1 * bmap

    # Safety check
    if len(bmap) < receptor_data_complete.shape[0]:
        raise ValueError(
            f"Brain map length ({len(bmap)}) is shorter than receptor map length "
            f"({receptor_data_complete.shape[0]}) in {dataset_name}"
        )

    # Use only first 218 rows to match receptor maps
    bmap = bmap[:218]

    # ------------------------------------------------------------------
    # Prepare results dictionary
    # ------------------------------------------------------------------
    results = {
        'receptor_set': [],
        'brain_map': [],
        'adj_r2': [],
        'moran_p': [],
        'regression_coef': []
    }

    # ------------------------------------------------------------------
    # Iterate through each receptor set
    # ------------------------------------------------------------------
    for set_name, receptors in receptor_sets.items():
        print(f"Analyzing receptor set: {set_name} -> {receptors}")

        idx = np.where(np.in1d(receptor_names, receptors))[0]
        receptor_data = receptor_data_complete[:, idx]

        if receptor_data.ndim == 1:
            receptor_data = receptor_data.reshape(-1, 1)

        # Regression
        r2, coef = get_reg_r_sq(zscore(receptor_data), zscore(bmap))
        print(f'    Adjusted R²: {r2:.4f}')
        print(f'    Coefficients: {coef}')

        # Spin-test
        surrogates = make_surrogates(bmap, dist_fname, spatnull, n_perm)
        pval = get_reg_r_pval(
            zscore(receptor_data),
            zscore(bmap),
            zscore(surrogates)
        )
        print(f'    Spin-test p-value: {pval:.4f}')

        # Store results
        results['receptor_set'].append(set_name)
        results['brain_map'].append(weight_col)
        results['adj_r2'].append(r2)
        results['moran_p'].append(pval)
        results['regression_coef'].append(
            coef.item() if np.size(coef) == 1 else coef.tolist()
        )

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# -------------------------------------------------------------------------
# Run analysis for joint CCA shared Mode 2
# -------------------------------------------------------------------------
datasets = [
    {
        "coef_file": project_dir + "results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/"
                                  "CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_coef.csv",
        "output_file": project_dir + "results/neurotransmitter/cmi/shared/"
                                    "individual_neurotransmitter_regression_results_shared_mode2.csv",
        "name": "CMI-shared-mode2"
    },
    {
        "coef_file": project_dir + "results/cca/stanford/wholebrain_cca_stanford_math_reading_combined/"
                                  "CCA_PCA_roi_gmv_brainnetome_mathreadingstd_ageinmodel_coef.csv",
        "output_file": project_dir + "results/neurotransmitter/stanford/shared/"
                                    "individual_neurotransmitter_regression_results_shared_mode2.csv",
        "name": "Stanford-shared-mode2"
    }
]

for ds in datasets:
    run_analysis_from_coef_csv(
        coef_file=ds["coef_file"],
        output_file=ds["output_file"],
        dataset_name=ds["name"],
        weight_col='xloading_m2',
        flip_sign=True
    )