# parseXNAT_Nazirah

This code will parse the BLSA data (updated: 2022) to calculate the statistics for each DTI measures in SLANT and EVE regions.

Atlases used:  
* EVE Type I, II, and III ([Oishi et al, 2009](https://pubmed.ncbi.nlm.nih.gov/19385016/)) 
* SLANT ([Huo et al, 2019](https://pubmed.ncbi.nlm.nih.gov/30910724/))

Software needed:  
* [MATLAB](https://www.mathworks.com/products/get-matlab.html)
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
* [MRtrix3](https://www.mrtrix.org)

## Running the parse code

Main code: MATLAB file `parseXNATdataWithADRD_nazirah_newdata.m`.  
* GitHub code path: https://github.com/MASILab/parseXNAT_Nazirah/tree/main/Nazirah
* MASI path: `/home/local/VANDERBILT/mohdkhn/Documents/Test3-newdata`  
(as of 2022/09/08)

Functions needed are in the `functions` folder.
<br/><br/>

Main thing needs to be changed is `D` and `D2`
* `D`: path to the data
* `D2`: path to codes (where the main code and functions directory are)

If you have different organization of the data, check all XXname variables if the paths are correct.
<br/><br/>

The original paths for each file used in this code is here:  
| File | Location |
|---|---|
|EVE Atlases | D2 > session_name > ASSESSORS > White_Matter_Stamper > WM_LABELS > EVEx_Labels.nii.gz |
|SLANT Atlas | D2 > session_name > ASSESSORS > Slant > REG > T1_seg.nii.gz |
|MPRAGE | D2 > session_name > SCANS > MPRAGE > NIFTI > XX_MPRAGE.nii.gz |
|FA/MD/RD/AD | D2 > session_name > ASSESSORS > DTI# > SCALARS > filename.nii.gz |

This code will generate images including (if not there yet):  
| File | Location |
|---|---|
|MPRAGE brain | session_name > SCANS > MPRAGE > NIFTI > mpr_brain.nii.gz
|MPRAGE brain mask | session_name > SCANS > MPRAGE > NIFTI > mpr_brain_mask.nii.gz
|XFM matrix from MPRAGE to B0 | session_name > SCANS > MPRAGE > NIFTI > mprage2b0.mat
|XFM matrix from B0 to MPRAGE (inv of above) | session_name > SCANS > MPRAGE > NIFTI > b02mprage.mat
|B0 image | D2 > session_name > ASSESSORS > DTI# > PREPROCESSED > b0.nii.gz
|all EVE2DTI images | D2 > session_name > ASSESSORS > White_Matter_Stamper > WM_LABELS > EVEx_Labels.nii.gz.subjLabels.nii.gz

<br/><br/>

## Additional info:   
[EVE and SLANT registration to FA/DWI](reg_method.md)  
[Working on ACCRE](work_on_accre.md)
