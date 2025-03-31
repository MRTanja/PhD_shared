#!/bin/bash

# This script processes DWI images with MRtrix3 to compute fixel-based metrics and perform whole-brain statistical analysis

# Input: preprocessed DWI images (after running dwidenoise, mrdegibbs for Gibbs ringing removal, 
# dwifslpreproc for motion and distortion correction, dwibiascorrect for bias field correction)

# Following steps in MRtrix3 documentation from here

#############################################################################################
#                                    GENERATING FODs                                        # 
############################################################################################# 

# Step 4: Computing (average) tissue response functions

base_dir="path/to/base/dir"
data_dir="$base_dir/data"
subset_dir="$base_dir/subset"
n= # number of threads to use 
design_matrix="/path/to/designmatrix.txt"
contrast_matrix="/path/to/contrastmatrix.txt"


cd data_dir
for_each * : dwi2response dhollander IN/*space-dwi_desc-preproc_dwi.mif IN/response_wm.txt IN/response_gm.txt IN/response_csf.txt

cd subset_dir
responsemean */response_wm.txt ../group_average_response_wm.txt
responsemean */response_gm.txt ../group_average_response_gm.txt
responsemean */response_csf.txt ../group_average_response_csf.txt

# Step 5: Upsampling DW images
cd data_dir
for_each -nthreads n * : mrgrid IN/*space-dwi_desc-preproc_dwi.mif regrid -vox 1.25 IN/dwi_preproc_upsampled.mif

# Step 6: Compute upsampled brain mask images
for_each -nthreads n * : dwi2mask IN/*dwi_preproc_upsampled.mif IN/dwi_mask_upsampled.mif

# Step 7: Fibre Orientation Distribution estimation (multi-tissue spherical deconvolution)
for_each -nthreads n * : dwi2fod msmt_csd IN/dwi_preproc_upsampled.mif ../group_average_response_wm.txt IN/wmfod.mif ../group_average_response_gm.txt IN/gm.mif  ../group_average_response_csf.txt IN/csf.mif -mask IN/dwi_mask_upsampled.mif

# Step 8: Joint bias field correction and intensity normalisation
for_each -nthreads n * : mtnormalise IN/wmfod.mif IN/wmfod_norm.mif IN/gm.mif IN/gm_norm.mif IN/csf.mif IN/csf_norm.mif -mask IN/dwi_mask_upsampled.mif


#############################################################################################
#                   GENERATING FOD TEMPLATE, REGISTERING ALL SUBJECTS                       # 
############################################################################################# 

# Step 9: Generate a study-specific unbiased FOD template
mkdir -p ../template/fod_input
mkdir ../template/mask_input

cd subset_dir
for_each * : ln -sr IN/wmfod_norm.mif ../template/fod_input/PRE.mif
for_each * : ln -sr IN/dwi_mask_upsampled.mif ../template/mask_input/PRE.mif
population_template ../template/fod_input -mask_dir ../template/mask_input ../template/wmfod_template.mif -voxel_size 1.25

# Step 10: Register all subject FOD images to the FOD template
cd data_dir
for_each -nthreads n * : mrregister IN/wmfod_norm.mif -mask1 IN/dwi_mask_upsampled.mif ../template/wmfod_template.mif -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif

# Step 11: Compute the template mask (intersection of all subject masks in template space)
for_each -nthreads n * : mrtransform IN/dwi_mask_upsampled.mif -warp IN/subject2template_warp.mif -interp nearest -datatype bit IN/dwi_mask_in_template_space.mif
mrmath */dwi_mask_in_template_space.mif min ../template/template_mask.mif -datatype bit

# Step 12: Compute a white matter template analysis fixel mask
fod2fixel -mask ../template/template_mask.mif -fmls_peak_value 0.06 ../template/wmfod_template.mif ../template/fixel_mask 

# Step 13: Warp FOD images to template space
for_each -nthreads n * : mrtransform IN/wmfod_norm.mif -warp IN/subject2template_warp.mif -reorient_fod no IN/fod_in_template_space_NOT_REORIENTED.mif

 
#############################################################################################
#                                    CALCULATE FD, LOGFC AND FDC                            # 
############################################################################################# 

# Step 14: Segment FOD images to estimate fixels and their apparent fibre density (FD)
for_each -nthreads n * : fod2fixel -mask ../template/template_mask.mif IN/fod_in_template_space_NOT_REORIENTED.mif IN/fixel_in_template_space_NOT_REORIENTED -afd fd.mif

# Step 15: Reorient fixels
for_each -nthreads n * : fixelreorient IN/fixel_in_template_space_NOT_REORIENTED IN/subject2template_warp.mif IN/fixel_in_template_space
for_each * : rm -r fixel_in_template_space_NOT_REORIENTED

# Step 16: Assign subject fixels to template fixels
for_each -nthreads n * : fixelcorrespondence IN/fixel_in_template_space/fd.mif ../template/fixel_mask ../template/fd PRE.mif

# Step 17: Compute the fibre cross-section (FC) metric
for_each -nthreads n * : warp2metric IN/subject2template_warp.mif -fc ../template/fixel_mask ../template/fc IN.mif

mkdir ../template/log_fc
cp ../template/fc/index.mif ../template/fc/directions.mif ../template/log_fc
for_each * : mrcalc ../template/fc/IN.mif -log ../template/log_fc/IN.mif

# Step 18: Compute a combined measure of fibre density and cross-section (FDC)
mkdir ../template/fdc
cp ../template/fc/index.mif ../template/fc/directions.mif ../template/fdc
for_each * : mrcalc ../template/fd/IN.mif ../template/fc/IN.mif -mult ../template/fdc/IN.mif


#############################################################################################
#                                    WHOLE-BRAIN ANALYSIS                                   # 
#############################################################################################

# Step 19: Perform whole-brain fibre tractography on the FOD template
cd ../template
tckgen -nthreads n -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image template_mask.mif -mask template_mask.mif -select 20000000 -cutoff 0.06 tracks_20_million.tck

# Step 20: Reduce biases in tractogram densities
tcksift tracks_20_million.tck wmfod_template.mif tracks_2_million_sift.tck -term_number 2000000

# Step 21: Generate fixel-fixel connectivity matrix
fixelconnectivity fixel_mask/ tracks_2_million_sift.tck matrix/

# Step 22: Smooth fixel data using fixel-fixel connectivity
fixelfilter fd smooth fd_smooth -matrix matrix/
fixelfilter log_fc smooth log_fc_smooth -matrix matrix/
fixelfilter fdc smooth fdc_smooth -matrix matrix/

# Step 23: Perform statistical analysis of FD, FC, and FDC
fixelcfestats fd_smooth/ files.txt design_matrix contrast_matrix matrix/ stats_fd/
fixelcfestats log_fc_smooth/ files.txt design_matrix contrast_matrix matrix/ stats_log_fc/
fixelcfestats fdc_smooth/ files.txt design_matrix contrast_matrix matrix/ stats_fdc/
