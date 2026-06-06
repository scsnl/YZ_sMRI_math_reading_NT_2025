# Glutamatergic signaling underlies brain structural organization for mathematical and reading abilities in children

© 2026 by The Board of Trustees of the Leland Stanford Junior University is licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)

For commercial license inquiries, contact the [Stanford Office of Technology Licensing (otl.stanford.edu)](https://otl.stanford.edu).

---

## Study Overview

The neurochemical mechanisms underlying individual differences in children's mathematical and reading abilities remain largely unknown. 
Here we investigate how neurotransmitter systems relate to brain structural organization supporting these abilities in two independent cohorts of children. 

Two key components:
1. **Brain Structural Phenotypes:** We derived whole-brain structural phenotypes associated with mathematical and reading abilities using structural MRI and Canonical Correlation Analysis (CCA).
2. **Neurotransmitter Mapping:** We mapped the spatial expression of these phenotypes onto the spatial distribution of neurotransmitter systems derived from PET data to identify which molecular systems are most strongly related to math- and reading-related brain structural organization.

We employed a **discovery-validation framework** using two independent cohorts (ages 7–14 years):
- **CMI-HBN discovery cohort (N = 760):** One of the largest neuroimaging studies of math and reading abilities.
- **Stanford validation cohort (N = 231):** Independent replication dataset.

By analyzing both **mathematical** and **reading** abilities in the same participants with identical methods, we directly compare molecular associations across cognitive domains.

---

## Project Structure

```
├── data/
│   ├── atlas/
│   │   ├── BN_Atlas_246_1mm.nii
│   │   ├── BN_Atlas_246_2mm.nii
│   │   ├── BNA_subregions.xlsx
│   │   ├── bn_atlas.xlsx
│   │   └── Shirer_14Networks/
│   ├── gmv/
│   │   ├── gmv_cmi_n760.csv
│   │   └── gmv_stanford_n231.csv
│   ├── subjectlist/
│   │   ├── subjectlist_cmi_n760.csv
│   │   └── subjectlist_stanford_n231.csv
│   └── neurotransmitter/
│       ├── PET_nifti_images/
│       ├── PET_parcellated/
│       │   └── bn246/
│       ├── nifti_images_19/
│       ├── receptor_data_bn246.csv
│       └── receptor_names_pet.npy
│
├── scripts/
│   ├── age_analysis/
│   ├── sex_analysis/
│   ├── behavior/
│   ├── BF/
│   ├── cca/
│   │   ├── full_sample/
│   │   ├── prediction/
│   │   ├── visualization/
│   │   └── table_top20.R
│   ├── control_analysis/
│   │   ├── cmi/
│   │   └── stanford/
│   │   ├── table_top20.R
│   │   ├── compare_gmv_weight_map_across_domains_cohorts.R
│   │   └── compare_gmv_weight_maps_across_models.R
│   ├── control_analysis/
│   │   ├── cmi/
│   │   ├── stanford/
│   │   └── partial_correlation.R
│   ├── network_analysis/
│   ├── neurotransmitter/
│   ├── shirer_to_bn/
│   ├── spin/
│   └── utility/
│
├── results/
│   ├── age_analysis/
│   ├── cca/
│   │   ├── cmi/
│   │   └── stanford/
│   ├── sex_analysis/
│   ├── cca/
│   │   ├── cmi/
│   │   ├── stanford/
│   │   ├── prediction/
│   │   ├── partial_corr/
│   │   └── brainmaps/
│   ├── network_analysis/
│   └── neurotransmitter/
│
└── README.md
```

**Note:** All scripts are thoroughly commented to guide users through each analysis step and make replication easier. After downloading the repository, users need to update the working directory in the R scripts (e.g., setwd) and the project directory in the Python scripts so that they point to the local path of the repository in the user’s environment. Once these repository-level paths are updated, the scripts, including interanll called scripts (e.g., myRfunc.R), should run smoothly.

---

## Input Data for Analyses

### Subjectlist
- **CMI-HBN subjectlist:** `subjectlist_cmi.csv` (N = 760)
- **Stanford subjectlist:** `subjectlist_stanford.csv` (N = 231)

### GMV (Gray Matter Volume)
- **CMI-HBN GMV:** `gmv_cmi_n760.csv` (rows = subjects, columns = brain regions)
- **Stanford GMV:** `gmv_stanford_n231.csv`

### Atlas
- **Brainnetome NIfTI atlas:** `BN_Atlas_246_2mm.nii`, `BN_Atlas_246_1mm.nii`
- **Brainnetome ROI details:** `bn_atlas.xlsx`
- **Brainnetome ROI MNI coordinates:** `BNA_subregions.xlsx`
- **Shirer 14 Networks:** NIfTI files under `Shirer_14Networks/`

### Neurotransmitter Data
- **Original PET images:** `data/neurotransmitter/PET_nifti_images/` (from [Hansen et al.](https://github.com/netneurolab/hansen_receptors/tree/main/data/PET_nifti_images))
- **Parcellated PET images:** `PET_parcellated/bn246/`
- **19 neurotransmitter maps (z-scored):** `nifti_images_19/`
- **Receptor matrix (region x neurotransmitter):** `receptor_data_bn246.csv`
- **Receptor names:** `receptor_names_pet.npy`

---

## Key Analysis Steps and Scripts

### 1. Behavioral Analysis
- **Script:** `scripts/behavior/behavior_analysis.R`

### 2. Canonical Correlation Analysis (CCA)
- **CMI-HBN:**  
  - Separate models: `scripts/cca/full_sample/Wholebrain_CCA_cmi.R`  
  - Seprate models' results: `results/cca/cmi/wholebrain_cca_cmi_math/`, `results/cca/cmi/wholebrain_cca_cmi_reading/`
  - Joint model: `scripts/cca/full_sample/Wholebrain_CCA_cmi_joint.R`  
  - Joint model's results: `results/cca/cmi/wholebrain_cca_cmi_math_reading_combined/`
  
- **Stanford:**  
  - Separate models: `scripts/cca/full_sample/Wholebrain_CCA_stanford.R`  
  - Seprate models' results: `results/cca/stanford/wholebrain_cca_stanford_math/`, `results/cca/stanford/wholebrain_cca_stanford_reading/`
  - Joint model: `scripts/cca/full_sample/Wholebrain_CCA_stanford_joint.R`  
  - Joint model's results: `results/cca/cmi/wholebrain_cca_stanford_math_reading_combined/`

- **Visualization:**  
  - Convert GMV weights to NIfTI: seperate models `scripts/cca/visualization/BN2brainmap_cca_SignChanged.m`, joint model `scripts/cca/visualization/BN2brainmap_cca_SignChanged_joint.m`
  - BrainNet Viewer for visualization.
- **Top 20% ROI Table:**  
  - `scripts/shirer_to_bn/map_shirer_to_brainnetome.m`  
  - `scripts/cca/table_top20.R`
- **Compare GMV weight maps:**  
  - across models: `scripts/cca/compare_gmv_weight_maps_across_models.m`  
  - across domains and cohorts: `scripts/cca/compare_gmv_weight_map_across_domains_cohorts.m`  
  
### 3. Prediction Analysis
- **Script:** `scripts/cca/prediction/prediction_cmi_model_to_stanford.R`
- **Scatter plot:** `scripts/cca/prediction/prediction_scatter_plot.R`

### 4. Age-Stratified Analysis
- **Script:** `scripts/age_analysis/ageAnalysis_with_plots.R`
- **Output:** GMV weights CSVs (full sample + 3 age bins).

### 5. Sex-Stratified Analysis
- **Script:** `scripts/sex_analysis/sexAnalysis_with_plots.R`
- **Output:** GMV weights CSVs (full sample + female/male).

### 5. Neurotransmitter Analysis
- **PET parcellation:** `scripts/neurotransmitter/parcellate.py`
- **Receptor matrix:** `scripts/neurotransmitter/make_receptor_matrix.py`
- **Spin test:** `scripts/spin/Calculating_Distancematrix_BN.m`
- **Regression:** 
  - separate models: `scripts/neurotransmitter/brainbeh_receptors_analysis.py`  
  - joint model: `scripts/neurotransmitter/brainbeh_receptors_analysis_joint.py`  
  (with `scripts/neurotransmitter/utility.py`)
- **FDR correction:** `scripts/add_fdrp.R`, `scripts/add_fdrp_joint.R`
- **Bayes Factor:**  
  - CMI-HBN: `scripts/BF/getBF_neurotransmitter.R`, `scripts/BF/getBF_neurotransmitter_joint.R`
  - Replication BF: `scripts/BF/getBFr_neurotransmitter.R`, `scripts/BF/getBFr_neurotransmitter_joint.R`
- **Visualization:** `scripts/neurotransmitter/simple_barplot_individual_receptors.R`, `scripts/neurotransmitter/simple_barplot_individual_receptors_joint.R`

### 6. Network Analysis
- `scripts/network_analysis/scatter_plot_directionality.R`
- `scripts/network_analysis/scatter_plot_directionality_network_shirer.R`

### 7. Control Analyses (IQ, SES, Site, TBV)
- **CMI-HBN:**  
  - IQ-controlled: `scripts/control_analysis/cmi/Wholebrain_CCA_cmi_ctrIQ.R`
  - SES-controlled: `scripts/control_analysis/cmi/Wholebrain_CCA_cmi_ctrSES.R`
  - Site-controlled: `scripts/control_analysis/cmi/Wholebrain_CCA_cmi_ctrSITE.R`
  - TVB-controlled: `scripts/control_analysis/partial_correlation.R`
  - Compare models: `scripts/control_analysis/cmi/check_weights_correlation.R`
- **Stanford:**  
  - IQ-controlled: `scripts/control_analysis/stanford/Wholebrain_CCA_stanford_ctrIQ.R`
  - SES-controlled: `scripts/control_analysis/stanford/controlAnalysis_SES_stanford.R`
  - TVB-controlled: `scripts/control_analysis/partial_correlation.R`
  - Compare models: `scripts/control_analysis/stanford/check_weights_correlation.R`
---

## How to Run
1. Clone this repository and ensure required R and Python dependencies are installed.
2. Update the working directory in the R scripts (e.g., setwd) and the project directory in the Python scripts
3. Prepare input data in `data/`.
4. Run behavioral and CCA scripts (e.g., `scripts/behavior`, `scripts/cca/full_sample`).
5. Use prediction, age-stratified, sex-stratified scripts for cross-validation.
6. Execute neurotransmitter and network analyses.
7. Visualize results using the provided R and MATLAB scripts.
8. All scripts include **detailed comments** to assist replication and adaptation.

---

