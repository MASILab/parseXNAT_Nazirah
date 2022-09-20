%% my code

%Files that we have:
D = '/Users/nana/Documents/MATLAB/BLSA/parseXNAT_Nazirah/Nazirah/testfileNFS';
subj = 'BLSA_5499';

D_subj = [D filesep subj];

%% Load label names
% T1 seg = same as andrew labels without CSF
T1seglabel = readcell([D_subj filesep 'T1_label_volumes.txt']);
T1SegLabelNames = T1seglabel(2:end,1);
T1SegLabelID = cell2mat(T1seglabel(2:end,2));


session = 'BLSA_5499_test';
D_session = [D_subj filesep session];

%% load DTI data
faname = [D_session filesep 'fa.nii.gz'];
mdname = [D_session filesep 'md.nii.gz'];
rdname = [D_session filesep 'rd.nii.gz'];
adname = [D_session filesep 'ad.nii.gz'];


fa = niftiread(faname);
md = niftiread(mdname);
rd = niftiread(rdname);
ad = niftiread(adname);

%% Find MPRAGE
mprfile = dir([D_session filesep '*MPRAGE.nii.gz']);
mprname = [D_session filesep mprfile(1).name];

%% Resampling
T1segname = [D_session filesep 'T1_seg.nii.gz'];

%Get tranform matrix
xfmname = [D_session 'mpr2fa.txt'];

if dir([D_session '*mpr2fa*'])<1
    system(['flirt -in ' mprname ' -out ' mprname '-flirt.nii.gz -ref ' faname ' -omat ' xfmname ' -dof 6'])
    system(['flirt -in ' T1segname ' -ref ' faname ' -applyxfm -init ' xfmname ' -interp nearestneighbour' ' -out ' label1name '-flirt.nii.gz'])
    system(['flirt -in ' T1segname ' -ref ' faname ' -applyxfm -init ' xfmname ' -interp nearestneighbour'])
end

% Resample the Multi-Atlas Labels: now how to get the transform file?
% T1seg is in MPRAGE size --> need to transform to DWI img size.
brainColorName = [T1segname '.subjLabels.nii.gz'];

% if resampled not exist yet, need to transform.
if(length(dir(brainColorName))<1)
    system(['reg_transform -ref ' faname ' -invAff ' xfmname ' ' xfmname '.inv']);
    system(['reg_resample -aff ' xfmname '.inv ' '-ref ' faname ' -flo ' T1segname ' -res ' T1segname '.subjLabels.nii.gz' ' -inter 0'])
end




