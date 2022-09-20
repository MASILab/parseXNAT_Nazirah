% clear all;

% Setup the root directory for all subjects
D = '/fs4/masi/landmaba/BLSAdti/BLSA/'

addpath('/fs4/masi/landmaba/BLSAdti/BLSA')
addpath(genpath('/fs4/masi/landmaba/BLSAdti/matlab'))

% Load the label names
EVElabelNames = textread([D 'EVE_Labels.csv'],'%s','delimiter','\n');
for i=1:length(EVElabelNames)
    EVElabelID(i) = sscanf(EVElabelNames{i},'%d');
end

BClabelNames = textread([D 'andrew_multiatlas_labels.csv'],'%s','delimiter','\n');
for i=1:length(BClabelNames)
    BClabelID(i) = sscanf(BClabelNames{i},'%d ',1);
end

DTISeglabelNames = textread([D 'DTIseg.txt'],'%s','delimiter','\n');
for i=1:length(DTISeglabelNames)
    DTISeglabelID(i) = sscanf(DTISeglabelNames{i},'%d ',1);
end
% Get a list of Subjects; ignore . and ..
SUBJS = dir([D filesep 'BLSA*']); SUBJS=SUBJS(1:end);

%% This is used to ensure that the stats file is replaced each time.
doMakePNG = 0;
doSkipPNGDone = 0;
doMakeReport = 1;

for jSubj = 1:length(SUBJS) %929 % [73 128 450 528 544 569 586 621 720 834 929]%834% 1:length(SUBJS)
    
    % Get a list of Sessions; ignore . and ..
    SESSIONS = dir([D SUBJS(jSubj).name]); SESSIONS =SESSIONS(3:end);
    
    for jSession=1:length(SESSIONS)
        
        DS = [D SUBJS(jSubj).name filesep SESSIONS(jSession).name filesep];
        Stamper = dir([DS  '*Stamper*']);
        
        ver3 =0;
        dtiQA = dir([DS '*dtiQA_v2_1']);
        if(length(dtiQA)~=2)
            dtiQA = dir([DS '*dtiQA_v2']);
        end
        if(length(dtiQA)~=2)
            dtiQA = dir([DS '*dtiQA_v2*']);
        end
        if((length(dtiQA)~=2)+strcmp(SESSIONS(jSession).name,'BLSA_4889_02-0_10')+strcmp(SESSIONS(jSession).name,'BLSA_4806_05-0_10'))
            ver3=1;
            dtiQA = dir([DS '*dtiQA_v3*']);
        end
        dtiQAMulti = dir([DS '*dtiQA_Multi*']);
        MultiAtlas = dir([DS '*Multi_Atlas*']);
        
        % Sometimes Multi-atlas will find another T1. If
        % more than one are found, choose the MPRAGE one.
        if(length(MultiAtlas)>1)
            MultiAtlas = dir([DS '*MPRAGE*Multi_Atlas*']);
        end
        MPRAGE = dir([DS 'MPRAGE*']);
        
        disp([SUBJS(jSubj).name ' ' SESSIONS(jSession).name])
        disp([length(Stamper) length(dtiQA) length(MultiAtlas)])
        try
            if(~and(and(length(dtiQA)>1,length(Stamper)==1),length(MultiAtlas)==1))
                % The data are not there. let's write a stats file that
                % tells the user why.
                if(and(and(length(dtiQA)<2,length(Stamper)==1),length(MultiAtlas)==1))
                    reportFileName = [D 'statsWithADRDVol' filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                    fp = fopen(reportFileName,'at');
                    fprintf(fp,'ONLY 1 DTI');
                    fclose(fp);
                    error('There are not 2 DTIs.');
                end
                
            else
                % if(and(and(length(dtiQA)>1,length(Stamper)==1),length(MultiAtlas)==1))
                if(length(dir([DS Stamper(1).name]))<3)
                    continue;
                end
                disp('FOUND valid dataset')
                
                reportFileName = [D 'statsWithADRDVol' filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                if(exist(reportFileName,'file'))
                    disp(['Already done: ' reportFileName]);
                    continue;
                end
                disp('touch')
                system(['echo `date` > ' reportFileName])
                
                if(doSkipPNGDone)
                    if(exist([D 'pngs' filesep SESSIONS(jSession).name '.png'],'file'))
                        disp('It''s done');
                        continue;
                    end
                end
                
                %% Find the MPRAGE
                mprfile = dir([DS MPRAGE(1).name filesep 'NIFTI' filesep '*.gz']);
                mprname = [DS MPRAGE(1).name filesep 'NIFTI' filesep mprfile(1).name];
                
                
                %% Deal with the Multi DTI session "DTI multi"
                cd([DS dtiQAMulti(1).name filesep 'TGZ'])
                files=dir();
                if(length(dir('QA_maps'))<3)
                    system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
                end
                
                adMname = [pwd filesep '..' filesep 'AD' filesep 'ad.nii.gz' ];
                rdMname = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                faMname= [pwd filesep 'QA_maps' filesep 'fa.nii'];
                mdMname= [pwd filesep 'QA_maps' filesep 'md.nii'];
                roiMname = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
                boxFABiasMname = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
                boxFAMname = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                boxFASigMname = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                
                verifyADRD(adMname,rdMname,faMname,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);
                
                
                %% Deal with the first DTI session "DTI(1)"
                if(ver3==0)
                    cd([DS dtiQA(1).name filesep 'TGZ'])
                    files=dir();
                    if(length(dir('QA_maps'))<3)
                        system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
                    end
                    
                    fa1name= [pwd filesep 'QA_maps' filesep 'fa.nii'];
                    md1name= [pwd filesep 'QA_maps' filesep 'md.nii'];
                    ad1name = [pwd filesep '..' filesep 'AD' filesep 'ad.nii.gz' ];
                    rd1name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                    roi1name = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
                    boxFABias1name = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
                    boxFA1name = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    boxFASig1name = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                    
                    verifyADRD(ad1name,rd1name,fa1name,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);
                else
                 
                    fa1name = [DS dtiQA(1).name filesep 'FA' filesep 'fa.nii.gz'];
                    
                    md1name= [DS dtiQA(1).name filesep 'MD' filesep 'md.nii.gz']; 
                    ad1name = [DS dtiQA(1).name filesep 'AD' filesep 'ad.nii.gz'];
                    rd1name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                    
                    roi1name = [DS dtiQA(1).name filesep 'extra' filesep  'multi_atlas_labels.nii.gz'];
                    boxFABias1name = NaN;%[DS dtiQA(1).name filesep 'extra' filesep 'BoxplotsBias.mat'];
                    boxFA1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    boxFASig1name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                end
                
                %% Resample the EVE and BrainColor masks from T1 to DTI(1) space
                masegname = [DS filesep MultiAtlas(1).name filesep 'SEG' filesep 'orig_target_seg.nii.gz'];
                
                % Use flirt to find a 6 dof transform
                %                 xfmname = [DS 'mpr2fa.txt'];
                %                 system(['flirt -in ' mprname ' -out ' mprname '-flirt.nii.gz -ref ' faMname ' -omat ' xfmname ' -dof 6'])
                %                 system(['flirt -in ' label1name ' -ref ' faMname ' -applyxfm -init ' xfmname ' -interp nearestneighbour' ' -out ' label1name '-flirt.nii.gz'])
                %                 system(['flirt -in ' label1name ' -ref ' faMname ' -applyxfm -init ' xfmname ' -interp nearestneighbour'])
                
                % Resample the FA labels from Eve into Subject FA
                xfm1name = [DS Stamper(1).name filesep '/TRANSFORMATION/fa_t1_transformation_matrix.txt'];
                if(length(dir(xfm1name))<1)
                    xfm1name = [DS Stamper(1).name filesep '/Intra_Session_Reg/outputAffine.txt'];
                end
                if(length(dir(xfm1name))<1)
                    error('Cannot find WM transform');
                end
                label1name = [DS Stamper(1).name filesep 'WM_LABELS' '/Rectified_EVE_Labels.nii.gz']
                eveName = [label1name '.subjLabels.nii.gz' ];
                if(length(dir(eveName))<1)
                    system(['reg_transform -ref ' fa1name ' -invAff ' xfm1name ' ' xfm1name '.inv']);
                    system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' label1name ' -res ' label1name '.subjLabels.nii.gz' ' -inter 0'])
                end
                
                
                % Resample the Multi-Atlas Labels
                brainColorName = [masegname '.subjLabels.nii.gz'];
                if(length(dir(brainColorName))<1)
                    system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' masegname ' -res ' masegname '.subjLabels.nii.gz' ' -inter 0'])
                end
                
                
                %% Now deal with the second DTI session
                if(ver3==0)
                    cd([DS dtiQA(2).name filesep 'TGZ'])
                    files=dir();
                    if(length(dir('QA_maps'))<3)
                        system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
                    end
                    ad2name = [pwd filesep '..' filesep 'AD' filesep 'ad.nii.gz' ];
                    rd2name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                    fa2name= [pwd filesep 'QA_maps' filesep 'fa.nii'];
                    md2name= [pwd filesep 'QA_maps' filesep 'md.nii'];
                    roi2name = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
                    boxFABias2name = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
                    boxFA2name = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    boxFASig2name = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                    
                    verifyADRD(ad2name,rd2name,fa2name,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);
                else
                    fa2name = [DS dtiQA(1).name filesep 'FA' filesep 'fa.nii.gz'];
                    
                    md2name= [DS dtiQA(1).name filesep 'MD' filesep 'md.nii.gz']; 
                    ad2name = [DS dtiQA(1).name filesep 'AD' filesep 'ad.nii.gz'];
                    rd2name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                    
                    roi2name = [DS dtiQA(1).name filesep 'extra' filesep  'multi_atlas_labels.nii.gz'];
                    boxFABias2name = NaN;%[DS dtiQA(1).name filesep 'extra' filesep 'BoxplotsBias.mat'];
                    boxFA2name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    boxFASig2name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                end
                %% Verify
                
                
                faM = loadniiorniigz([faMname]);
                fa1 = loadniiorniigz([fa1name ]);
                fa2 = loadniiorniigz([fa2name ]);
                mdM = loadniiorniigz([mdMname ]);
                md1 = loadniiorniigz([md1name ]);
                md2 = loadniiorniigz([md2name ]);
                
                adM = loadniiorniigz(adMname);
                ad1 = loadniiorniigz(ad1name);
                ad2 = loadniiorniigz(ad2name);
                rdM = loadniiorniigz(rdMname);
                rd1 = loadniiorniigz(rd1name);
                rd2 = loadniiorniigz(rd2name);
                eve = loadniiorniigz(eveName);
                bc = loadniiorniigz(brainColorName);
                
                roi1 = loadniiorniigz([roi1name ]);
                roi2 = loadniiorniigz([roi2name ]);
                
                
                % debug problem with differing nifti headers (fixed in DTI qa 2.1)
                if(mean(mean(mean((fa1.img-faM.img).^2)))>mean(mean(mean((flipdim(fa1.img,2)-faM.img).^2))))
                    disp('Flipping');
                    fa1.img = flipdim(fa1.img,2);
                    md1.img = flipdim(md1.img,2);
                    ad1.img = flipdim(ad1.img,2);
                    rd1.img = flipdim(rd1.img,2);
                    
                    fp=fopen([D 'FatalErrorsWithADRD.txt'],'at');
                    fprintf(fp,'WARN: %s - %s\n',SESSIONS(jSession).name,'flipping dti1');
                    fclose(fp);
                end
                
                % debug problem with differing nifti headers (fixed in DTI qa 2.1)
                if(mean(mean(mean((fa2.img-faM.img).^2)))>mean(mean(mean((flipdim(fa2.img,2)-faM.img).^2))))
                    disp('Flipping');
                    fa2.img = flipdim(fa2.img,2);
                    md2.img = flipdim(md2.img,2);
                    ad2.img = flipdim(ad2.img,2);
                    rd2.img = flipdim(rd2.img,2);
                    fp=fopen([D 'FatalErrorsWithADRD.txt'],'at');
                    fprintf(fp,'WARN: %s - %s\n',SESSIONS(jSession).name,'flipping dti2');
                    fclose(fp);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Figures
                
                if(doMakeReport)
                    
                    ColHeader = {'Session'};
                    ColValues = {SESSIONS(jSession).name};
                    
                    %% PG 4 - Load stats (bias, variance of FA and MD by ROI's)
                    %                 boxFAMname = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    %                 boxFABiasMname = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
                    
                    %                 boxFASigMname = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                    
                    try
                    boxFA1 = load(boxFA1name);
                    boxFA2 = load(boxFA2name);
                    for j=1:length(DTISeglabelNames)
                        ColHeader{end+1} = ['FAMED-DTI1-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = nanmedian(boxFA1.FAroi(boxFA1.grp==(2*j-1)));
                        ColHeader{end+1} = ['FAMED-DTI2-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = nanmedian(boxFA2.FAroi(boxFA2.grp==(2*j-1)));
                    end
                    catch
                       
                        for j=1:length(DTISeglabelNames)
                        ColHeader{end+1} = ['FAMED-DTI1-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = NaN; 
                        ColHeader{end+1} = ['FAMED-DTI2-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = NaN; 
                        end
                    
                    end
                        
                    
                    try
                    boxFA1 = load(boxFASig1name);
                    boxFA2 = load(boxFASig2name);
                    for j=1:length(DTISeglabelNames)
                        ColHeader{end+1} = ['FASIG-DTI1-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = nanmedian(boxFA1.FAbootROI(boxFA1.grp==(2*j-1)));
                        ColHeader{end+1} = ['FASIG-DTI2-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = nanmedian(boxFA2.FAbootROI(boxFA2.grp==(2*j-1)));
                    end
                    catch
                         for j=1:length(DTISeglabelNames)
                        ColHeader{end+1} = ['FASIG-DTI1-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = NaN;
                        ColHeader{end+1} = ['FASIG-DTI2-ROI-' DTISeglabelNames{j}];
                        ColValues{end+1} = NaN;
                    end
                    end
                    
                    try
                        boxFA1 = load(boxFABias1name);
                        boxFA2 = load(boxFABias2name);
                        for j=1:length(DTISeglabelNames)
                            ColHeader{end+1} = ['FABIAS-DTI1-ROI-' DTISeglabelNames{j}];
                            ColValues{end+1} = nanmedian(boxFA1.Biasroi(boxFA1.grp==(2*j-1)));
                            ColHeader{end+1} = ['FABIAS-DTI2-ROI-' DTISeglabelNames{j}];
                            ColValues{end+1} = nanmedian(boxFA2.Biasroi(boxFA2.grp==(2*j-1)));
                        end
                    catch
                        for j=1:length(DTISeglabelNames)
                            ColHeader{end+1} = ['FABIAS-DTI1-ROI-' DTISeglabelNames{j}];
                            ColValues{end+1} = NaN;
                            ColHeader{end+1} = ['FABIAS-DTI2-ROI-' DTISeglabelNames{j}];
                            ColValues{end+1} = NaN;
                        end
                    end
                    
                    %% PG 6 - DSC on labels
                    for j=1:length(DTISeglabelNames)
                        ColHeader{end+1} = ['STATS-ROI-' DTISeglabelNames{j} '-DSC'];
                        A = roi1.img==j;
                        B = roi2.img==j;
                        ColValues{end+1} = 2 * sum(A(:).*B(:)) / (sum(A(:))+sum(B(:)));
                    end
                    
                    %% PG 8, 10, 12, 14 FA, MD by Eve and BrainColor label
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-FA-mean'];
                        ColValues{end+1} = mean(fa1.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-FA-std'];
                        ColValues{end+1} = std(fa1.img(eve.img(:)==EVElabelID(j)));
                        
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-FA-mean'];
                        ColValues{end+1} = mean(fa2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-FA-std'];
                        ColValues{end+1} = std(fa2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-FA-mean'];
                        ColValues{end+1} = mean(faM.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-FA-std'];
                        ColValues{end+1} = std(faM.img(eve.img(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-MD-mean'];
                        ColValues{end+1} = mean(md1.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-MD-std'];
                        ColValues{end+1} = std(md1.img(eve.img(:)==EVElabelID(j)));
                        
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-MD-mean'];
                        ColValues{end+1} = mean(md2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-MD-std'];
                        ColValues{end+1} = std(md2.img(eve.img(:)==EVElabelID(j)));
                        
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-MD-mean'];
                        ColValues{end+1} = mean(mdM.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-MD-std'];
                        ColValues{end+1} = std(mdM.img(eve.img(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-FA-mean'];
                        ColValues{end+1} = mean(fa1.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-FA-std'];
                        ColValues{end+1} = std(fa1.img(bc.img(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-FA-mean'];
                        ColValues{end+1} = mean(fa2.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-FA-std'];
                        ColValues{end+1} = std(fa2.img(bc.img(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-FA-mean'];
                        ColValues{end+1} = mean(faM.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-FA-std'];
                        ColValues{end+1} = std(faM.img(bc.img(:)==BClabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-MD-mean'];
                        ColValues{end+1} = mean(md1.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-MD-std'];
                        ColValues{end+1} = std(md1.img(bc.img(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-MD-mean'];
                        ColValues{end+1} = mean(md2.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-MD-std'];
                        ColValues{end+1} = std(md2.img(bc.img(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-MD-mean'];
                        ColValues{end+1} = mean(mdM.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-MD-std'];
                        ColValues{end+1} = std(mdM.img(bc.img(:)==BClabelID(j)));
                    end
                    %%%%%%%%%%%%%%% AD
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-AD-mean'];
                        ColValues{end+1} = mean(ad1.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-AD-std'];
                        ColValues{end+1} = std(ad1.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-AD-mean'];
                        ColValues{end+1} = mean(ad2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-AD-std'];
                        ColValues{end+1} = std(ad2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-AD-mean'];
                        ColValues{end+1} = mean(adM.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-AD-std'];
                        ColValues{end+1} = std(adM.img(eve.img(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-AD-mean'];
                        ColValues{end+1} = mean(ad1.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-AD-std'];
                        ColValues{end+1} = std(ad1.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-AD-mean'];
                        ColValues{end+1} = mean(ad2.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-AD-std'];
                        ColValues{end+1} = std(ad2.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-AD-mean'];
                        ColValues{end+1} = mean(adM.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-AD-std'];
                        ColValues{end+1} = std(adM.img(bc.img(:)==BClabelID(j)));
                    end
                    
                    
                    %%%%%%%%%%%%%%% RD
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-RD-mean'];
                        ColValues{end+1} = mean(rd1.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-RD-std'];
                        ColValues{end+1} = std(rd1.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-RD-mean'];
                        ColValues{end+1} = mean(rd2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-RD-std'];
                        ColValues{end+1} = std(rd2.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-RD-mean'];
                        ColValues{end+1} = mean(rdM.img(eve.img(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-RD-std'];
                        ColValues{end+1} = std(rdM.img(eve.img(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-RD-mean'];
                        ColValues{end+1} = mean(rd1.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-RD-std'];
                        ColValues{end+1} = std(rd1.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-RD-mean'];
                        ColValues{end+1} = mean(rd2.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-RD-std'];
                        ColValues{end+1} = std(rd2.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-RD-mean'];
                        ColValues{end+1} = mean(rdM.img(bc.img(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-RD-std'];
                        ColValues{end+1} = std(rdM.img(bc.img(:)==BClabelID(j)));
                    end
                    
                    %%%%%%%%%%%%%%% ROI Volumes
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'Volume'];
                        ColValues{end+1} = sum(bc.img(:)==BClabelID(j))*prod(bc.hdr.dime.pixdim(2:4));
                    end
                    
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'Volume'];
                        ColValues{end+1} = sum(eve.img(:)==EVElabelID(j))*prod(eve.hdr.dime.pixdim(2:4));
                    end
                    
                    %%%%%%%%%%%%%%%  Write out
                    
                    if(~exist([D 'AllStats-HeaderWithADRDVol.csv'],'file'))
                        fp = fopen([D 'AllStats-HeaderWithADRDVol.csv'],'wt'); fprintf(fp,'%s, ',ColHeader{:}); fprintf(fp,'\n'); fclose(fp)
                    end
                    fp = fopen(reportFileName,'at'); fprintf(fp,'%s, ',ColValues{1}); fprintf(fp,'%e, ',ColValues{2:end}); fprintf(fp,'\n'); fclose(fp)
                end
                %% Make a preview figure
                if(doMakePNG)
                    r=find(fa1.img(:,128,32));
                    x = round([min(r) median(r) max(r)]);
                    fa1.img(x,:,:)=1;
                    faM.img(x,:,:)=1;
                    eve.img(x,:,:)=128;
                    bc.img(x,:,:)=128;
                    figure(1);
                    clf;
                    imagesc(flipud([eve.img(:,:,32) faM.img(:,:,32)*255; bc.img(:,:,32) fa1.img(:,:,32)*255; fa2.img(:,:,32)*255 abs(fa1.img(:,:,32)-fa2.img(:,:,32))*255]') ); axis equal tight off
                    title([SUBJS(jSubj).name ' ' SESSIONS(jSession).name])
                    
                    drawnow
                    saveas(gcf,[D 'pngs' filesep SESSIONS(jSession).name '.png'])
                end
                
            end
        catch err
            disp('Oh no!')
            fp=fopen([D 'FatalErrorsWithADRDv7.txt'],'at');
            fprintf(fp,'%s - %s\n',SESSIONS(jSession).name,err.message);
            fclose(fp);
        end
        
    end
end
