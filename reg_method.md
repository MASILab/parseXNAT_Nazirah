EVE and SLANT registration to FA/DWI

0. **Get b0**

`dwiextract $dwi - -bzero -fslgrad $bvec $bval -quiet | mrmath - mean $b0 -axis 3 -quiet`
1. **Brain extract the MPRAGE**

`bet $mprage $mpragebrain -R -f 0.5 -g 0 -m`

2. **DTI to MPRAGE brain (mpragebrain) --> will get b02mprage.mat**

`epi_reg --epi=$b0 --t1=$mprage --t1brain=$mpragebrain --out=$b02mprage`

3. **Convert [DTI to MPRAGE brain] to [MPRAGE to DTI]**

`convert_xfm -omat mprage2b0.mat -inverse b02mprage.mat`

4. **EVE (label1name) to FA (fa1name) using mprage2b0.mat**

`flirt -in $label1name -ref $fa1name -applyxfm -init mprage2b0.mat -out eve2fa.nii.gz -interp nearestneighbour`

5) **SLANT (slant) to FA (fa1name) usng mprage2b0.mat**

`flirt -in $slant -ref $fa1name -applyxfm -init mprage2b0.mat -out slant2fa.nii.gz -interp nearestneighbour`

**Software used**
* FSL
* MRtrix3

**NOTE:**
* Be sure to check on every image produced by these commands to be what you are expecting
* Best to register using b0 to register to mpragebrain
* How to get b0? use bval and bvec to find where bval=0 and take the mean of the ones in the dwmri.nii.gz
