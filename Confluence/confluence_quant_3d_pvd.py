#################################################################################################
#                   Confluence quantification 3D for deep and periventricular WM                #
#                                   Tanja Schmidt, 25.05.24                                     #
#################################################################################################

# This script takes WMH segmentations (binarized or non-binarized) and calculates a 
# confluence metric for each subject, in 3D, separately for periventricular WM and deep WM

# Changes you'll need to make so that script works with your paths/filenames: lines 63, 69 and 73

import numpy as np
import pandas as pd
import nibabel as nib
import glob

# Function calculate_confluence calculates the metric for one slice
s = 0.05 # Defines width of Gaussian kernel, this was the optimal value in my tests

def calculate_confluence(image,s):
    # Get the x, y and z coordinates of each voxel in the image
    x, y, z = np.indices(image.shape)
    x_coord = x.flatten()
    y_coord = y.flatten()
    z_coord = z.flatten()
    vox_values = image.flatten()
    coord = pd.DataFrame({
        'vox_value': vox_values,
        'x': x_coord,
        'y': y_coord,
        'z': z_coord
    })
    # Delete all 0 voxels
    coord = coord[coord['vox_value']!=0].reset_index()
    # Create df with products of elementwise multiplication (multiply every voxel with every voxel)
    products = np.outer(coord['vox_value'],coord['vox_value'])
    products = pd.DataFrame(products)
    # Elementwise multiplication yields redundant values (a*b = b*a, a*a) -> set those to 0
    products[:] = np.tril(products.values, k=-1)
    equ_input = products.stack().reset_index()
    equ_input.columns = ['row', 'col', 'product']
    # Get rid of values that were set to 0
    equ_input = equ_input[equ_input['product']!=0]
    # Identify i, j, m and n coordinates that need to enter equation in next step
    equ_input['i'] = coord.loc[equ_input['row'],'x'].to_numpy()
    equ_input['j'] = coord.loc[equ_input['row'],'y'].to_numpy()
    equ_input['k'] = coord.loc[equ_input['row'],'z'].to_numpy()
    equ_input['m'] = coord.loc[equ_input['col'],'x'].to_numpy()
    equ_input['n'] = coord.loc[equ_input['col'],'y'].to_numpy()
    equ_input['o'] = coord.loc[equ_input['col'],'z'].to_numpy()
    # Calculate all elements of confluence metric
    equ_input['result'] = equ_input['product'] * np.exp(-s*(abs(equ_input['i']-equ_input['m'])**2 + abs(equ_input['j']-equ_input['n'])**2 + abs(equ_input['k']-equ_input['o'])**2))
    # Sum all elements to calculate confluence metric for the volume
    confluence = equ_input['result'].sum()
    return(confluence)

# Function calculate_volume calculates the number of WMH voxels in a volume
def calculate_volume(image):
    volume = np.count_nonzero(image)
    return volume


# Load images, run function calculate_confluence while looping through subjects and slices:
base_dir = '/home/ts887/rds/hpc-work/BIANCA/BIANCA_output/' # Change to your directory that contains images
WM = ['d','pv'] # Deep and periventricular white matter
input_dict = {}
for wm in WM:
    confluence_df = pd.DataFrame()
    volume_df = pd.DataFrame()
    for sub in glob.glob(base_dir + f'*thr06_{wm}.nii.gz'): # Change '*thr06*' to string that all image filenames contain
        print(f'Processing subject {sub}')
        image = nib.load(sub).get_fdata()
        # Get subject ID from filename -> adapt this so it works with your file names, replace '9_' and '_thr' with strings left and right of subject ID
        sub_id = sub.split('9_', 1)[1].split('_thr',1)[0]
        confluence_val = calculate_confluence(image, s)
        volume_val = calculate_volume(image)
        confluence_df = pd.concat([confluence_df, pd.DataFrame({'WBIC_ID': [f'{sub_id}'], f'confluence_{wm}': [confluence_val]})], ignore_index=True)
        volume_df = pd.concat([volume_df, pd.DataFrame({'WBIC_ID': [f'{sub_id}'], f'volume_{wm}': [volume_val]})], ignore_index=True)
  
    # Output so far: confluence_df and volume_df = dataframes with one row per subject,
    # containing the confluence metric for that slice and the number of WMH voxels
    confluence_final = pd.merge(confluence_df, volume_df, on='WBIC_ID', how='inner')
    # Normalize confluence metric with WMH volume
    confluence_final[f'confluence_norm_{wm}'] = confluence_final[f'confluence_{wm}']/confluence_final[f'volume_{wm}']
    # Normalize with maximum possible confluence value, i.e. value for an image where every voxel is a WMH (240.5 is the value for matrix size 192*256*256)
    confluence_final[f'confluence_scaled_{wm}'] = confluence_final[f'confluence_norm_{wm}']/240.5
    input_dict[wm] = confluence_final

result = pd.merge(input_dict['d'],input_dict['pv'],on=['WBIC_ID'])

# result = dataframe with one row per subject, contains all metrics (I left everything in for sanity checks),
# periventricular WM and deep WM in one dataframe,
# the relevant one is in the last column: confluence_norm_scaled, this will be a value between 0 and 1
