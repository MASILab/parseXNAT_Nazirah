if(ver3==0)
    [ad1name,rd1name,fa1name,md1name,roi1name,boxFABias1name,boxFA1name,boxFASig1name] = get_and_verify_ADRD([DS dtiQA(1).name filesep 'TGZ']);
else
    
    fa1name = [DS dtiQA(1).name filesep 'FA' filesep 'fa.nii'];
    md1name= [DS dtiQA(1).name filesep 'MD' filesep 'md.nii'];
    ad1name = [DS dtiQA(1).name filesep 'AD' filesep 'ad.nii'];
    rd1name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii' ];
    
    roi1name = [DS dtiQA(1).name filesep 'extra' filesep  'multi_atlas_labels.nii.gz'];
    boxFABias1name = NaN;%[DS dtiQA(1).name filesep 'extra' filesep 'BoxplotsBias.mat'];
    boxFA1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
    boxFASig1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
end