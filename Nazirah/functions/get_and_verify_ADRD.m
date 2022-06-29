function [adname,rdname,faname,mdname,roiname,boxFABiasname,boxFAname,boxFASigname] = get_and_verify_ADRD(location)

cd(location)
files=dir();
if(length(dir('QA_maps'))<3)
    system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
end

% Get path for AD,RD,FA,MD, etc files
adname = [pwd filesep '..' filesep 'AD' filesep 'ad1.nii.gz' ];
rdname = [pwd filesep '..' filesep 'RD' filesep 'rd1.nii.gz' ];
faname= [pwd filesep 'QA_maps' filesep 'fa.nii'];
mdname= [pwd filesep 'QA_maps' filesep 'md.nii'];
roiname = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
boxFABiasname = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
boxFAname = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
boxFASigname = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];

% Create all AD, RD folder
verifyADRD2(adname,rdname,faname,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);