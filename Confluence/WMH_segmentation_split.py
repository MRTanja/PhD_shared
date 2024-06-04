#################################################################################################
#                        Separation of WMH into deep and periventricular                        #
#                                   Tanja Schmidt, 25.05.24                                     #
#################################################################################################


# This script takes WMH segmentations and splits each into 2 segmentations: one with 
# periventricular WMH and one with deep WMH


# From FSL BIANCA userguide (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BIANCA/Userguide#BiancaMasking):

#################################################################################################################################################

# To extract WMH volumes within a certain distance from the ventricles, as a measure of periventricular (or deep) WMH, you can do the following:

    # 1. Generate a ventricle mask (<structural image>_vent.nii.gz) from your T1 image with make_bianca_mask (see section Masking for details)
    # (if needed) register it to your base image

    # 2. Create a distance map from the ventricles using FSL function distancemap, where the intensity of each voxel represent its distance in mm 
    # from the ventricles mask: distancemap <structural image>_vent.nii.gz dist_to_vent

    # 3. Threshold the distance map to obtain the region within <dist> mm of the ventricles, fslmaths dist_to_vent -uthr <dist> -bin dist_to_vent_periventricular or, 
    # alternatively, to look for deep WMH (more than <dist> mm from the ventricles) fslmaths dist_to_vent -thr <dist> -bin dist_to_vent_deep

#################################################################################################################################################

# In order to obtain <structural image>_vent.nii.gz: apply make_bianca_mask to T1 images, and for that you will first
# need to run fsl_anat on T1 images to create required input: 
# CSF partial volume estimation (*pve_0.nii.gz) and warp file MNI2structural (*MNI_to_T1_nonlin_field.nii.gz)

#################################################################################################################################################


import os
import glob

base_dir = '/home/ts887/rds/hpc-work/BIANCA/BIANCA_images/CUH_BIDS/'


# Run fsl_anat on T1 images to create input necessary for make_bianca_mask
for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/'):
    t1 = glob.glob(sub + '*rec-norm_T1w.nii.gz')[0].split('.nii.gz', 1)[0]
    com = 'fsl_anat --nosubcortseg -i ' + t1
    print(com)
    res = os.system(com)


# Create ventmask with make_bianca_mask and input from fsl_anat
for anat in glob.glob(base_dir + 'sub-*/ses-*/anat/*T1w.anat/'):
    t1_biascorr = glob.glob(anat + '*biascorr.nii.gz')[0].split('.nii.gz', 1)[0]
    pve = glob.glob(anat + '*pve_0.nii.gz')[0].split('.nii.gz', 1)[0]
    mni_to_t1_field = glob.glob(anat + '*MNI_to_T1_nonlin_field.nii.gz')[0].split('.nii.gz', 1)[0]
    com = 'make_bianca_mask ' + t1_biascorr + ' ' + pve + ' ' + mni_to_t1_field + ' 0'
    print(com)
    res = os.system(com)


# Register ventmask to FLAIR
for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/'):
    flair_candidate = glob.glob(sub + '*rec-norm_FLAIR.nii.gz')
    t1_candidate = glob.glob(sub + '*anat/T1.nii.gz')
    if flair_candidate and t1_candidate: # Check if all necessary files exist
        flair = flair_candidate[0].split('.nii.gz', 1)[0]
        t1 = t1_candidate[0].split('.nii.gz', 1)[0]
        com = 'flirt -in ' + t1 + ' -ref ' + flair + ' -out ' + t1 + '_to_FLAIR -omat ' + t1 + '_to_FLAIR.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp nearestneighbour'
        print(com)
        res = os.system(com)
    else:
        print(f'No matching file found in {sub}')


for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/'):
    flair_candidate = glob.glob(sub + '*rec-norm_FLAIR.nii.gz')
    ventmask_candidate = glob.glob(sub + '*anat/*ventmask.nii.gz')
    mat_candidate = glob.glob(sub + '*anat/T1_to_FLAIR.mat')
    print(mat_candidate)
    if flair_candidate and ventmask_candidate and mat_candidate:
        flair = flair_candidate[0].split('.nii.gz', 1)[0]
        ventmask = ventmask_candidate[0].split('.nii.gz', 1)[0]
        mat = mat_candidate[0]
        com = 'flirt -in ' + ventmask + '-applyxfm -init ' + mat + ' -out ' + ventmask + '_to_FLAIR ' + '-paddingsize 0.0 -interp nearestneighbour -ref ' + flair
        print(com)
        res = os.system(com)
    else:
        print(f'No matching file found in {sub}')


# Run distancemap on registered ventricle masks, output: *dist_to_vent.nii.gz where intensity of each voxel represent its distance in mm from ventricles mask
for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/*T1w.anat/'):
    ventmask_candidate = glob.glob(sub + '*ventmask_to_FLAIR.nii.gz')
    if ventmask_candidate:
        ventmask = ventmask_candidate[0].split('.nii.gz', 1)[0]
        com = 'distancemap -i ' + ventmask + ' -o ' + ventmask + '_dist_to_vent'
        print(com)
        res = os.system(com)
    else:
        print(f'No matching file found in {sub}')


# Threshold distance map to create two separate maps (periventricular WM = within 5 voxels around ventricles, deep WM = outside of that)
for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/*T1w.anat/'):
    distmap_candidate = glob.glob(sub + '*dist_to_vent.nii.gz')
    if distmap_candidate:
        distmap = distmap_candidate[0].split('.nii.gz', 1)[0]
        com = 'fslmaths ' + distmap + ' -uthr 5 -bin ' + distmap + '_pv'
        print(com)
        res = os.system(com)
        com = 'fslmaths ' + distmap + ' -thr 5 -bin ' + distmap + '_d'
        print(com)
        res = os.system(com)
    else:
        print(f'No matching file found in {sub}')


# Add periventricular WM mask to ventricle mask to create new pv WM mask (to make sure no WMH voxels at the very edge of ventricles are left out)
for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/*T1w.anat/'):
    pv_candidate = glob.glob(sub + '*_vent_pv.nii.gz')
    ventmask_candidate = glob.glob(sub + '*ventmask_to_FLAIR.nii.gz')
    if pv_candidate and ventmask_candidate:
        pv = pv_candidate[0].split('.nii.gz', 1)[0]
        ventmask = ventmask_candidate[0].split('.nii.gz',1)[0]
        com = 'fslmaths ' + pv + ' -add ' + ventmask + ' ' + pv + '_full'
        print(com)
        res = os.system(com)
    else:
        print(f'No matching file found in {sub}')


# Apply deep WM mask (d) and periventricular WM masks (pv) to BIANCA segmentation, or to the WMH segmentation that you have
for sub in glob.glob(base_dir + 'sub-*/ses-*/anat/'):
    id = sub.split('sub-')[-1].split('/ses')[0]
    bianca_candidate = glob.glob(f'/home/ts887/rds/hpc-work/BIANCA/BIANCA_CUH_output/*{id}_thr06.nii.gz') # Use WMH segmentation here
    pv_candidate = glob.glob(sub + '*anat/*vent_pv_full.nii.gz')
    d_candidate = glob.glob(sub + '*anat/*vent_d.nii.gz')
    if bianca_candidate and pv_candidate and d_candidate:
        bianca = bianca_candidate[0].split('.nii.gz', 1)[0]
        pv = pv_candidate[0].split('.nii.gz', 1)[0]
        d = d_candidate.[0].split('.nii.gz', 1)[0]
        com = 'fslmaths ' + bianca + ' -mas ' + pv + ' ' + bianca + '_pv'
        print(com)
        res = os.system(com)
        com = 'fslmaths ' + bianca + ' -mas ' + d + ' ' + bianca + '_d'
        print(com)
        res = os.system(com)
    else:
        print(f'No matching file found in {sub}')


# Images ending in *_pv.nii.gz and *_d.nii.gz are input on which we'll calculate confluence metric separately for periventricular and deep WM

