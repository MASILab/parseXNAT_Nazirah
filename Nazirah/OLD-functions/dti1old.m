if(ver3==0)
    cd([DS dtiQA(1).name filesep 'TGZ'])
    files=dir();
    if(length(dir('QA_maps'))<3)
        system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
    end
    
    fa1name= [pwd filesep 'QA_maps' filesep 'fa.nii'];
    md1name= [pwd filesep 'QA_maps' filesep 'md.nii'];
    ad1name = [pwd filesep '..' filesep 'AD' filesep 'ad2.nii.gz' ];
    rd1name = [pwd filesep '..' filesep 'RD' filesep 'rd2.nii.gz' ];
    roi1name = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
    boxFABias1name = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
    boxFA1name = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
    boxFASig1name = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
    
    verifyADRD(ad1name,rd1name,fa1name,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);
else
    
    fa1name = [DS dtiQA(1).name filesep 'FA' filesep 'fa.nii.gz'];
    
    md1name= [DS dtiQA(1).name filesep 'MD' filesep 'md.nii.gz'];
    ad1name = [DS dtiQA(1).name filesep 'AD' filesep 'ad2.nii.gz'];
    rd1name = [pwd filesep '..' filesep 'RD' filesep 'rd2.nii.gz' ];
    
    roi1name = [DS dtiQA(1).name filesep 'extra' filesep  'multi_atlas_labels.nii.gz'];
    boxFABias1name = NaN;%[DS dtiQA(1).name filesep 'extra' filesep 'BoxplotsBias.mat'];
    boxFA1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
    boxFASig1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
end