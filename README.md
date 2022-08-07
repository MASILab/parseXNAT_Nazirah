# parseXNAT_Nazirah

This code will parse the BLSA data (updated: 2022) to calculate the statistics for each DTI measures in SLANT and EVE regions.

Software needed:
* [MATLAB](https://www.mathworks.com/products/get-matlab.html)
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
* [MRtrix3](https://www.mrtrix.org)

Main code: MATLAB file `parseXNATdataWithADRD_nazirah_newdata.m`.

Functions needed are in the `functions` folder.

Main thing needs to be changed is `D` and `D2`, which are the path to the data and the code.
If you have different organization of the data, check all XXname variables if the paths are correct.

The original paths for each file used in this code is here:
* Atlases: D > atlas_name
* MPRAGE: D2 > session_name > SCANS > MPRAGE > NIFTI > XX_MPRAGE.nii.gz
* FA/MD/RD/AD: D2 > session_name > ASSESSORS > DTI# > SCALARS > filename.nii.gz

This code will generate few files (if not there yet):
* MPRAGE brain: session_name > SCANS > MPRAGE > NIFTI > mpr_brain.nii.gz
* MPRAGE brain mask: session_name > SCANS > MPRAGE > NIFTI > mpr_brain_mask.nii.gz
* XFM matrix from MPRAGE to B0: session_name > SCANS > MPRAGE > NIFTI > mprage2b0.mat
* XFM matrix from B0 to MPRAGE (inv of above): session_name > SCANS > MPRAGE > NIFTI > b02mprage.mat
* B0 image: D2 > session_name > ASSESSORS > DTI# > PREPROCESSED > b0.nii.gz
* Registered SLANT2FA image: D2 > session_name > ASSESSORS > slant > SEG > T1_seg.nii.gz.subjLabels.nii.gz
* Registered EVE2FA image: D2 > session_name > ASSESSORS > White_Matter_Stamper > WM_LABELS > Rectified_EVE_Labels.nii.gz.subjLabels.nii.gz


Processes (step-by-step in the main code)

To be listed here

[EVE and SLANT registration to FA/DWI](reg_method.md)
[Working on ACCRE](work_on_accre.md)
