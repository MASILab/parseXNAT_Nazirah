cd([DS dtiQAMulti(1).name filesep 'TGZ'])
files=dir();
if(length(dir('QA_maps'))<3)
    system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
end

% Get path for AD,RD,FA,MD, etc files
adMname = [pwd filesep '..' filesep 'AD' filesep 'ad.nii.gz' ];
rdMname = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
faMname= [pwd filesep 'QA_maps' filesep 'fa.nii'];
mdMname= [pwd filesep 'QA_maps' filesep 'md.nii'];
roiMname = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
boxFABiasMname = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
boxFAMname = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
boxFASigMname = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];

% Create all AD, RD folder
verifyADRD2(adMname,rdMname,faMname,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);