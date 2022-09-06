%ONLY FOR TESTING
% 2022/06/23: cannot test since BLSA_5499 already have ad.nii.gz

%cd([DS dtiQAMulti(1).name filesep 'TGZ'])
files=dir();
if(length(dir('QA_maps'))<3)
    system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
end

adMname = [pwd filesep '..' filesep 'AD' filesep 'ad2.nii.gz' ];
rdMname = [pwd filesep '..' filesep 'RD' filesep 'rd2.nii.gz' ];
faMname= [pwd filesep 'QA_maps' filesep 'fa.nii'];
mdMname= [pwd filesep 'QA_maps' filesep 'md.nii'];
roiMname = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
boxFABiasMname = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
boxFAMname = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
boxFASigMname = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];

verifyADRD(adMname,rdMname,faMname,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);