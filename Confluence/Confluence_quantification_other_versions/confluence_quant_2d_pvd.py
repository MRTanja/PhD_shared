#################################################################################################
#                   Confluence quantification for deep and periventricular WM                   #
#                                   Tanja Schmidt, 14.02.24                                     #
#################################################################################################

# This script takes WMH segmentations (binarized or non-binarized) and calculates a 
# confluence metric for each subject, separately for deep and periventricular WM

# Changes you'll need to make so that script works with your paths/filenames: lines 58, 66 and 70

import numpy as np
import pandas as pd
import nibabel as nib
import glob
from functools import reduce

# Function calculate_confluence calculates the metric for one slice
s = 0.05 # Defines width of Gaussian kernel, this was the optimal value in my tests

def calculate_confluence(image,slice,s):
    # Pick one 2D slice from 3D numpy array, convert into pandas dataframe:
    slice_data = pd.DataFrame(image[:,:,slice]) 
    # Use only slices with non-zero voxels, no need to process slices with no WMH voxels:
    if slice_data.any(axis=None): 
        # Create df with voxel coordinates and voxel values
        coord = slice_data.stack().reset_index()
        coord.columns =['x', 'y', 'vox_value'] 
        # Delete all 0 voxels
        coord = coord[coord['vox_value']!=0].reset_index() 
        # Create df with products of elementwise multiplication (multiply every voxel with every voxel)
        coord_vox_value = coord['vox_value'].to_numpy()
        products = pd.DataFrame(coord_vox_value * coord_vox_value[:,None])
        # Elementwise multiplication yields redundant values (a*b = b*a, a*a) -> set those to 0
        products[:] = np.tril(products.values, k=-1)
        equ_input = products.stack().reset_index()
        equ_input.columns = ['row', 'col', 'product']
        # Get rid of values that were set to 0
        equ_input = equ_input[equ_input['product']!=0]
        # Identify i, j, m and n coordinates that need to enter equation in next step
        equ_input['i'] = coord.loc[equ_input['row'],'x'].to_numpy()
        equ_input['j'] = coord.loc[equ_input['row'],'y'].to_numpy()
        equ_input['m'] = coord.loc[equ_input['col'],'x'].to_numpy()
        equ_input['n'] = coord.loc[equ_input['col'],'y'].to_numpy()
        # Calculate all elements of confluence metric
        equ_input['result'] = equ_input['product'] * np.exp(-s*(abs(equ_input['i']-equ_input['m'])**2 + abs(equ_input['j']-equ_input['n'])**2))
        # Sum all elements to calculate confluence metric for the slice
        confluence = equ_input['result'].sum()
        return(confluence)

# Function calculate_volume calculates the number of WMH voxels in a slice
def calculate_volume(image,slice):
    slice_data = pd.DataFrame(image[:,:,slice])
    volume = slice_data.sum().sum()
    return volume


# Load images, run function calculate_confluence while looping through subjects and slices:
base_dir = '/home/ts887/rds/hpc-work/BIANCA/BIANCA_output/' # Change to your directory that contains images

WM = ['d','pv'] # Deep and periventricular white matter; 
# this assumes you have two sets of WMH segmentations, one ending in *_d.nii.gz for hyperintensities in deep WM 
# and one ending in *_pv.nii.gz for periventricular WM; the script WMH_segmentation_split.py generates these
# two sets of segmentations from normal whole-brain WMH segmentations

input_dict = {}

for wm in WM:
    confluence_df = pd.DataFrame()
    volume_df = pd.DataFrame()
    for sub in glob.glob(base_dir + f'*thr06_{wm}.nii.gz'): # Change '*thr06*' to string that all image filenames contain
        print(f'Processing subject {sub}')
        image = nib.load(sub).get_fdata() 
        # Get subject ID from filename -> adapt this so it works with your file names, replace '9_' and '_thr' with strings left and right of subject ID
        sub_id = sub.split('9_', 1)[1].split('_thr',1)[0]
        confluence_list = [calculate_confluence(image, slice, s) for slice in range(image.shape[2])]
        volume_list = [calculate_volume(image, slice) for slice in range(image.shape[2])]
        confluence_df[f'{sub_id}'] = pd.DataFrame(confluence_list)
        volume_df[f'{sub_id}'] = pd.DataFrame(volume_list)

    # For each slice, normalize confluece metric by dividing it by number of WMH voxels in that slice
    conf_norm = confluence_df/volume_df

    # Calculate sums across all slices for each subject, create dataframe
    conf_sum = pd.DataFrame({'WBIC_ID': confluence_df.sum().index, f'confluence_{wm}': confluence_df.sum().values})
    vol_sum = pd.DataFrame({'WBIC_ID': volume_df.sum().index, f'volume_{wm}': volume_df.sum().values})
    conf_norm_sum = pd.DataFrame({'WBIC_ID': conf_norm.sum().index, f'confluence_norm_{wm}': conf_norm.sum().values})
    slice_count = confluence_df.apply(lambda x: (x > 0).sum())
    nonzero_slices = pd.DataFrame({'WBIC_ID': confluence_df.columns, f'nonzero_slices_{wm}': slice_count}).reset_index(drop=True)

    dfs = [conf_sum, vol_sum, conf_norm_sum, nonzero_slices]
    confluence_final = reduce(lambda left,right: pd.merge(left,right,on=['WBIC_ID'],
                                                how='inner'), dfs)
    confluence_final[f'confluence_scaled_{wm}'] = confluence_final[f'confluence_norm_{wm}']/(30.20349728*confluence_final[f'nonzero_slices_{wm}'])
    confluence_final['WBIC_ID'] = confluence_final['WBIC_ID'].astype(int)

    input_dict[wm] = confluence_final

result = pd.merge(input_dict['d'],input_dict['pv'],on=['WBIC_ID'])
# result = dataframe with one row per subject, contains all metrics (I left everything in for sanity checks),
# the relevant one is in the last column: Confluence_norm_scaled, this will be a value between 0 and 1
