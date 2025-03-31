# -*- coding: utf-8 -*-
"""
Parcellate volumetric PET images
"""

import numpy as np
# from nilearn.datasets import fetch_atlas_schaefer_2018
from neuromaps.parcellate import Parcellater

scale = 'bn246'

brainnetome = '/Users/zhangyuan/Google Drive/2020_LongitudinalStructure_Reading_Math/BN_Atlas_246_2mm.nii.gz'

path = "/Users/zhangyuan/Google Drive/2020_LongitudinalStructure_Reading_Math/hansen_receptors-main/data/PET_nifti_images/"
outpath = "/Users/zhangyuan/Google Drive/2020_LongitudinalStructure_Reading_Math/hansen_receptors-main/data/PET_parcellated/"+scale+"/"

receptors_nii = [path+'5HT1a_way_hc36_savli.nii',
                 path+'5HT1a_cumi_hc8_beliveau.nii',
                 path+'5HT1b_az_hc36_beliveau.nii',
                 path+'5HT1b_p943_hc22_savli.nii',
                 path+'5HT1b_p943_hc65_gallezot.nii.gz',
                 path+'5HT2a_cimbi_hc29_beliveau.nii',
                 path+'5HT2a_alt_hc19_savli.nii',
                 path+'5HT2a_mdl_hc3_talbot.nii.gz',
                 path+'5HT4_sb20_hc59_beliveau.nii',
                 path+'5HT6_gsk_hc30_radhakrishnan.nii.gz',
                 path+'5HTT_dasb_hc100_beliveau.nii',
                 path+'5HTT_dasb_hc30_savli.nii',
                 path+'A4B2_flubatine_hc30_hillmer.nii.gz',
                 path+'CB1_omar_hc77_normandin.nii.gz',
                 path+'CB1_FMPEPd2_hc22_laurikainen.nii',
                 path+'D1_SCH23390_hc13_kaller.nii',
                 path+'D2_fallypride_hc49_jaworska.nii',
                 path+'D2_flb457_hc37_smith.nii.gz',
                 path+'D2_flb457_hc55_sandiego.nii.gz',
                 path+'D2_raclopride_hc7_alakurtti.nii',
                 path+'DAT_fpcit_hc174_dukart_spect.nii',
                 path+'DAT_fepe2i_hc6_sasaki.nii.gz',
                 path+'GABAa-bz_flumazenil_hc16_norgaard.nii',
                 path+'GABAa_flumazenil_hc6_dukart.nii',
                 path+'H3_cban_hc8_gallezot.nii.gz',
                 path+'M1_lsn_hc24_naganawa.nii.gz',
                 path+'mGluR5_abp_hc22_rosaneto.nii',
                 path+'mGluR5_abp_hc28_dubois.nii',
                 path+'mGluR5_abp_hc73_smart.nii',
                 path+'MU_carfentanil_hc204_kantonen.nii',
                 path+'MU_carfentanil_hc39_turtonen.nii',
                 path+'NAT_MRB_hc77_ding.nii.gz',
                 path+'NAT_MRB_hc10_hesse.nii',
                 path+'NMDA_ge179_hc29_galovic.nii.gz',
                 # path+'VAChT_feobv_hc3_spreng.nii',
                 path+'VAChT_feobv_hc4_tuominen.nii',
                 path+'VAChT_feobv_hc5_bedard_sum.nii',
                 path+'VAChT_feobv_hc18_aghourian_sum.nii']

parcellated = {}
parcellater = Parcellater(brainnetome, 'MNI152')

for receptor in receptors_nii:
    parcellated[receptor] = parcellater.fit_transform(receptor, 'MNI152', True)
    name = receptor.split('/')[-1]  # get nifti file name
    name = name.split('.')[0]  # remove .nii
    np.savetxt(outpath+name+'.csv', parcellated[receptor], delimiter=',')
