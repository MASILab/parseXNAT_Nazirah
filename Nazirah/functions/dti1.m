if(ver3==0)
    [ad1name,rd1name,fa1name,md1name,roi1name,boxFABias1name,boxFA1name,boxFASig1name] = get_and_verify_ADRD([DS dtiQA(1).name filesep 'TGZ']);
else
    
    fa1name = [DS dtiQA(1).name filesep 'FA' filesep 'fa.nii.gz'];
    
    md1name= [DS dtiQA(1).name filesep 'MD' filesep 'md.nii.gz'];
    ad1name = [DS dtiQA(1).name filesep 'AD' filesep 'ad1.nii.gz'];
    rd1name = [pwd filesep '..' filesep 'RD' filesep 'rd1.nii.gz' ];
    
    roi1name = [DS dtiQA(1).name filesep 'extra' filesep  'multi_atlas_labels.nii.gz'];
    boxFABias1name = NaN;%[DS dtiQA(1).name filesep 'extra' filesep 'BoxplotsBias.mat'];
    boxFA1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
    boxFASig1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
end