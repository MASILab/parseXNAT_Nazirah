% Setup the root directory for all subjects

%% Original PATH
%D = '/fs4/masi/landmaba/BLSAdti/BLSA/'
%addpath('/fs4/masi/landmaba/BLSAdti/BLSA')
%addpath(genpath('/fs4/masi/landmaba/BLSAdti/matlab'))

%% Nazirah's Mac
% D = '/Users/nana/Documents/MATLAB/BLSA/original_files/';
% D2 = '/Users/nana/Documents/MATLAB/BLSA/parseXNAT_Nazirah/';
% 
% addpath([D2 'Nazirah/functions/']);

%% MASI Computers on Nazirah's Mac
D = '/Users/nana/masi-42/Documents/data/';
D2 = '/Users/nana/masi-42/Documents/Test1/';
addpath([D2 'functions/']);

%% MASI Computers
% D = '/home/local/VANDERBILT/mohdkhn/Documents/data/'; % data path
% D2 = '/home/local/VANDERBILT/mohdkhn/Documents/Test1/'; % code and labels path
% addpath([D2 'functions/']);


%% set time to record date and time to add to filename (FatalErrors.txt, AllStats.csv)
timedate = datestr(now,'yyyymmddTHHMMSS');
ReportFolder = [D 'statsWithADRDVol' timedate];

if ~exist(ReportFolder,'dir')
    mkdir(ReportFolder)
end

FatalErrorFile = [D 'FatalErrorsWithADRDv8-' timedate '.txt'];
AllStatsFile = [D 'AllStats-HeaderWithADRDVol-' timedate '.csv'];

%%
% Load the label names
% EVE_path = [D 'EVE_Labels.csv'];
% BC_path = [D 'andrew_multiatlas_labels.csv'];
% DTI_path = [D 'DTIseg.txt'];

EVE_path = [D2 'EVE_Labels.csv'];
BC_path = [D2 'andrew_multiatlas_labels.csv'];
DTI_path = [D2 'DTIseg.txt'];

[EVElabelNames,EVElabelID,BClabelNames,BClabelID,DTISeglabelNames,DTISeglabelID] = get_label_names(EVE_path,BC_path,DTI_path);

%% Get a list of Subjects; ignore . and ..
%SUBJS = dir([D filesep 'BLSA*']); SUBJS=SUBJS(1:end);
SUBJS = dir([D 'BLSA*']);


%% This is used to ensure that the stats file is replaced each time.
doMakePNG = 0;
doSkipPNGDone = 0;
doMakeReport = 1;

for jSubj = 1:length(SUBJS) %929 % [73 128 450 528 544 569 586 621 720 834 929]%834% 1:length(SUBJS)
    
    % Get a list of Sessions; ignore . and ..
    %SESSIONS = dir([D SUBJS(jSubj).name]); SESSIONS =SESSIONS(3:end);
    SESSIONS = dir([D SUBJS(jSubj).name]);
    SESSIONS(startsWith({SESSIONS.name},'.')) = [];
    
    fprintf('%s\n', SUBJS(jSubj).name)
    if ~strcmp(SUBJS(jSubj).name,'BLSA_4889')
        continue;
    end
    
    for jSession=1:length(SESSIONS)
        DS = [D SUBJS(jSubj).name filesep SESSIONS(jSession).name filesep];
        Stamper = dir([DS  '*Stamper*']);
        
        ver3 = 0;
        dtiQA = dir([DS '*dtiQA_v2_1']);
        if(length(dtiQA)~=2)
            dtiQA = dir([DS '*dtiQA_v2']);
            fprintf('Got dtiQA_v2\n')
        end
        
        if(length(dtiQA)~=2)
            dtiQA = dir([DS '*dtiQA_v2*']);
            fprintf('Got *dtiQA_v2*\n')
        end
        
        if((length(dtiQA)~=2)+strcmp(SESSIONS(jSession).name,'BLSA_4889_02-0_10')+strcmp(SESSIONS(jSession).name,'BLSA_4806_05-0_10'))
            ver3 = 1;
            dtiQA = dir([DS '*dtiQA_v3*']);
            fprintf('Got *dtiQA_v3*\n')
        end
        dtiQAMulti = dir([DS '*dtiQA_Multi*']);
        MultiAtlas = dir([DS '*Multi_Atlas*']);
        
        % Sometimes Multi-atlas will find another T1. If
        % more than one are found, choose the MPRAGE one.
        if(length(MultiAtlas)>1)
            MultiAtlas = dir([DS '*MPRAGE*Multi_Atlas*']);
            fprintf('More than 1 multiatlas\n')
        end
        MPRAGE = dir([DS 'MPRAGE*']);
        
        fprintf([SUBJS(jSubj).name ' ' SESSIONS(jSession).name '\n']);
        fprintf('Stamper: %d dtiQA: %d MultiAtlas: %d\n', length(Stamper),length(dtiQA), length(MultiAtlas));
        
        try
            % RULE:
            % dtiQA has to be more than 1
            % WM stamper must exist (equals to 1)
            % MultiAtlas must exist (equals to 1)
            
            % IF dtiQA < 2 and HAS Stamper AND MultiAtlas
            %if(~and(and(length(dtiQA)>1,length(Stamper)==1),length(MultiAtlas)==1))
            if ~(length(dtiQA)>1 && length(Stamper)==1 && length(MultiAtlas)==1)
                % The data are not there. let's write a stats file that
                % tells the user why.
                
                % if(and(and(length(dtiQA)<2,length(Stamper)==1),length(MultiAtlas)==1))
                if length(dtiQA)<2 && length(Stamper)==1 && length(MultiAtlas)==1
                    % This write report only if dtiQA is less than 2.
                    % Stamper and MultiAtlas must exist.
                    reportFileName = [ReportFolder filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                    fp = fopen(reportFileName,'at');
                    fprintf(fp,'ONLY 1 DTI');
                    fclose(fp);
                    error('There are not 2 DTIs.');
                end
                
                % Now, we have enough DTI files, let's check Stamper, CSV, and PNG files
            else
                % if(and(and(length(dtiQA)>1,length(Stamper)==1),length(MultiAtlas)==1))
                
                % CHECK: if WM Stamper folder has less than 3 files, don't
                % continue. Move to the next session.
                if(length(dir([DS Stamper(1).name]))<3)
                    continue;
                end
                fprintf('FOUND valid dataset\n')
                
                % Create CSV report file.
                reportFileName = [ReportFolder filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                
                % If CSV file exist, don't continue below. Move to the next Session
                if(exist(reportFileName,'file'))
                    fprintf(['Already done: ' reportFileName '\n']);
                    continue;
                end
                fprintf('touch\n')
                system(['echo `date` > ' reportFileName]);
                
                % If PNG exists, don't continue below and go to next session. [Whyy?]
                if(doSkipPNGDone)
                    if(exist([D 'pngs' filesep SESSIONS(jSession).name '.png'],'file'))
                        fprintf('It''s done\n');
                        continue;
                    end
                end
                
                %% NOW, ALL FILES EXISTS AND NEED TO BE RUN
                %% Find the MPRAGE
                mprfile = dir([DS MPRAGE(1).name filesep 'NIFTI' filesep '*.gz']);
                mprname = [DS MPRAGE(1).name filesep 'NIFTI' filesep mprfile(1).name];
                
                
                %% Deal with the Multi DTI session "DTI multi"
                % In TGZ folder...
                %multiDTI
                [adMname,rdMname,faMname,mdMname,roiMname,boxFABiasMname,boxFAMname,boxFASigMname] = get_and_verify_ADRD([DS dtiQAMulti(1).name filesep 'TGZ']);
                
                
                %% Deal with the first DTI session "DTI(1)"
                dti1
                % if ver3 == 1, don't run verifyADRD
                
                %% Resample the EVE and BrainColor masks from T1 to DTI(1) space
                masegname = [DS MultiAtlas(1).name filesep 'SEG' filesep 'orig_target_seg.nii.gz'];
                
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
                label1name = [DS Stamper(1).name filesep 'WM_LABELS' '/Rectified_EVE_Labels.nii.gz'];
                eveName = [label1name '.subjLabels.nii.gz' ];
                if(length(dir(eveName))<1)
                    
                    % reg_transform -ref ref_img -invAff transform_mat transform_mat.inv
                    % link: http://cmictig.cs.ucl.ac.uk/wiki/index.php/Reg_transform
                    system(['reg_transform -ref ' fa1name ' -invAffine ' xfm1name ' ' xfm1name '.inv']);
                    % input: xfm1name
                    % output: xfm1name.inv
                    %system(['convert_xfm -omat ' xfm1name '2.inv' ' -inverse ' xfm1name]);
                    
                    % reg_resample -aff transform_mat.inv -ref ref_img -flo EVE_Labels.nii.gz -res EVE_Labels.nii.gz.subjLabels.nii.gz -inter 0
                    % link: http://cmictig.cs.ucl.ac.uk/wiki/index.php/Reg_resample
                    % input: xfm1name, label1name, fa1name
                    % output: label1name.subjLabels.nii.gz
                    system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' label1name ' -res ' label1name '.subjLabels.nii.gz' ' -inter 0'])
                    %system(['flirt -in ' label1name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' label1name '.subjLabels.nii.gz'])
                end
                
                
                % Resample the Multi-Atlas Labels
                brainColorName = [masegname '.subjLabels.nii.gz'];
                if(length(dir(brainColorName))<1)
                    system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' masegname ' -res ' masegname '.subjLabels.nii.gz' ' -inter 0'])
                end
                
                
                %% Now deal with the second DTI session
                if(ver3==0)
                    %                     cd([DS dtiQA(2).name filesep 'TGZ'])
                    %                     files=dir();
                    %                     if(length(dir('QA_maps'))<3)
                    %                         system(['tar xvf ' files(3).name ' --exclude=*Reg*'])
                    %                     end
                    %                     ad2name = [pwd filesep '..' filesep 'AD' filesep 'ad.nii.gz' ];
                    %                     rd2name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                    %                     fa2name= [pwd filesep 'QA_maps' filesep 'fa.nii.gz'];
                    %                     md2name= [pwd filesep 'QA_maps' filesep 'md.nii.gz'];
                    %                     roi2name = [pwd filesep 'extra' filesep 'multi_atlas_labels.nii'];
                    %                     boxFABias2name = [pwd filesep 'extra' filesep 'BoxplotsBias.mat'];
                    %                     boxFA2name = [pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    %                     boxFASig2name = [pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                    %
                    %                     verifyADRD2(ad2name,rd2name,fa2name,[pwd filesep 'QA_maps' filesep 'dt.Bdouble']);
                    
                    [ad2name,rd2name,fa2name,md2name,roi2name,boxFABias2name,boxFA2name,boxFASig2name] = get_and_verify_ADRD([DS dtiQA(2).name filesep 'TGZ']);
                else
                    %                     fa2name = [DS dtiQA(1).name filesep 'FA' filesep 'fa.nii.gz'];
                    %                     md2name = [DS dtiQA(1).name filesep 'MD' filesep 'md.nii.gz'];
                    %                     ad2name = [DS dtiQA(1).name filesep 'AD' filesep 'ad.nii.gz'];
                    %                     rd2name = [pwd filesep '..' filesep 'RD' filesep 'rd.nii.gz' ];
                    %
                    %                     roi2name = [DS dtiQA(1).name filesep 'extra' filesep  'multi_atlas_labels.nii.gz'];
                    
                    ad2name = findfileniiorgz([DS dtiQA(1).name filesep 'AD'], 'ad.nii');
                    rd2name = findfileniiorgz([pwd filesep '..' filesep 'RD'], 'rd.nii');
                    fa2name = findfileniiorgz([DS dtiQA(1).name filesep 'FA'],'fa.nii');
                    md2name = findfileniiorgz([DS dtiQA(1).name filesep 'MD'],'md.nii');
                    roi2name = findfileniiorgz([DS dtiQA(1).name filesep 'extra'], 'multi_atlas_labels.nii');
                    
                    boxFABias2name = NaN;%[DS dtiQA(1).name filesep 'extra' filesep 'BoxplotsBias.mat'];
                    boxFA2name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFA.mat'];
                    boxFASig2name = NaN;%[pwd filesep 'extra' filesep 'BoxplotsFAsigma.mat'];
                end
                %% Verify
                
                %faM = loaduntouchniiorniigz(faMname);
                %faM2 = niftiread(faMname);
                
                faM = loadniiorgz(faMname); faMinfo = infoniiorgz(faMname);
                fa1 = loadniiorgz(fa1name); fa1info = infoniiorgz(fa1name);
                fa2 = loadniiorgz(fa2name); fa2info = infoniiorgz(fa2name);
                
                mdM = loadniiorgz(mdMname); %mdMinfo = niftiinfo(mdMname);
                md1 = loadniiorgz(md1name); %md1info = niftiinfo(md1name);
                md2 = loadniiorgz(md2name); %md2info = niftiinfo(md2name);
                
                adM = loadniiorgz(adMname);
                ad1 = loadniiorgz(ad1name);
                ad2 = loadniiorgz(ad2name);
                
                rdM = loadniiorgz(rdMname);
                rd1 = loadniiorgz(rd1name);
                rd2 = loadniiorgz(rd2name);
                
                eve = loadniiorgz(eveName); eveinfo = infoniiorgz(eveName);
                bc = loadniiorgz(brainColorName); bcinfo = infoniiorgz(brainColorName);
                
                roi1 = loadniiorgz(roi1name);
                roi2 = loadniiorgz(roi2name);
                
                % debug problem with differing nifti headers (fixed in DTI qa 2.1)
                if(mean(mean(mean((fa1-faM).^2)))>mean(mean(mean((flip(fa1,2)-faM).^2))))
                    fprintf('Flipping\n');
                    fa1 = flip(fa1,2);
                    md1 = flip(md1,2);
                    ad1 = flip(ad1,2);
                    rd1 = flip(rd1,2);
                    
                    fp=fopen(FatalErrorFile,'at');
                    fprintf(fp,'WARN: %s - %s\n',SESSIONS(jSession).name,'flipping dti1');
                    fclose(fp);
                end
                
                % debug problem with differing nifti headers (fixed in DTI qa 2.1)
                if(mean(mean(mean((fa2-faM).^2)))>mean(mean(mean((flip(fa2,2)-faM).^2))))
                    fprintf('Flipping\n');
                    fa2 = flip(fa2,2);
                    md2 = flip(md2,2);
                    ad2 = flip(ad2,2);
                    rd2 = flip(rd2,2);
                    fp=fopen(FatalErrorFile,'at');
                    fprintf(fp,'WARN: %s - %s\n',SESSIONS(jSession).name,'flipping dti2');
                    fclose(fp);
                end
                
                
                fp = fopen([D 'headerCheckv8.txt'],'at');
                fprintf(fp,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',SESSIONS(jSession).name,...
                    bcinfo.PixelDimensions(1),bcinfo.PixelDimensions(2),bcinfo.PixelDimensions(3),...
                    eveinfo.PixelDimensions(1),eveinfo.PixelDimensions(2),eveinfo.PixelDimensions(3),...
                    fa1info.PixelDimensions(1),fa1info.PixelDimensions(2),fa1info.PixelDimensions(3),...
                    fa2info.PixelDimensions(1),fa2info.PixelDimensions(2),fa2info.PixelDimensions(3),...
                    faMinfo.PixelDimensions(1),faMinfo.PixelDimensions(2),faMinfo.PixelDimensions(3));
                fclose(fp);
                
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
                        A = roi1==j;
                        B = roi2==j;
                        ColValues{end+1} = 2 * sum(A(:).*B(:)) / (sum(A(:))+sum(B(:)));
                    end
                    
                    %% PG 8, 10, 12, 14 FA, MD by Eve and BrainColor label
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-FA-mean'];
                        ColValues{end+1} = mean(fa1(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-FA-std'];
                        ColValues{end+1} = std(fa1(eve(:)==EVElabelID(j)));
                        
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-FA-mean'];
                        ColValues{end+1} = mean(fa2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-FA-std'];
                        ColValues{end+1} = std(fa2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-FA-mean'];
                        ColValues{end+1} = mean(faM(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-FA-std'];
                        ColValues{end+1} = std(faM(eve(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-MD-mean'];
                        ColValues{end+1} = mean(md1(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-MD-std'];
                        ColValues{end+1} = std(md1(eve(:)==EVElabelID(j)));
                        
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-MD-mean'];
                        ColValues{end+1} = mean(md2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-MD-std'];
                        ColValues{end+1} = std(md2(eve(:)==EVElabelID(j)));
                        
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-MD-mean'];
                        ColValues{end+1} = mean(mdM(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-MD-std'];
                        ColValues{end+1} = std(mdM(eve(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-FA-mean'];
                        ColValues{end+1} = mean(fa1(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-FA-std'];
                        ColValues{end+1} = std(fa1(bc(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-FA-mean'];
                        ColValues{end+1} = mean(fa2(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-FA-std'];
                        ColValues{end+1} = std(fa2(bc(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-FA-mean'];
                        ColValues{end+1} = mean(faM(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-FA-std'];
                        ColValues{end+1} = std(faM(bc(:)==BClabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-MD-mean'];
                        ColValues{end+1} = mean(md1(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-MD-std'];
                        ColValues{end+1} = std(md1(bc(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-MD-mean'];
                        ColValues{end+1} = mean(md2(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-MD-std'];
                        ColValues{end+1} = std(md2(bc(:)==BClabelID(j)));
                        
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-MD-mean'];
                        ColValues{end+1} = mean(mdM(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-MD-std'];
                        ColValues{end+1} = std(mdM(bc(:)==BClabelID(j)));
                    end
                    %%%%%%%%%%%%%%% AD
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-AD-mean'];
                        ColValues{end+1} = mean(ad1(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-AD-std'];
                        ColValues{end+1} = std(ad1(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-AD-mean'];
                        ColValues{end+1} = mean(ad2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-AD-std'];
                        ColValues{end+1} = std(ad2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-AD-mean'];
                        ColValues{end+1} = mean(adM(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-AD-std'];
                        ColValues{end+1} = std(adM(eve(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-AD-mean'];
                        ColValues{end+1} = mean(ad1(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-AD-std'];
                        ColValues{end+1} = std(ad1(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-AD-mean'];
                        ColValues{end+1} = mean(ad2(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-AD-std'];
                        ColValues{end+1} = std(ad2(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-AD-mean'];
                        ColValues{end+1} = mean(adM(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-AD-std'];
                        ColValues{end+1} = std(adM(bc(:)==BClabelID(j)));
                    end
                    
                    
                    %%%%%%%%%%%%%%% RD
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-RD-mean'];
                        ColValues{end+1} = mean(rd1(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI1-RD-std'];
                        ColValues{end+1} = std(rd1(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-RD-mean'];
                        ColValues{end+1} = mean(rd2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTI2-RD-std'];
                        ColValues{end+1} = std(rd2(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-RD-mean'];
                        ColValues{end+1} = mean(rdM(eve(:)==EVElabelID(j)));
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'DTIM-RD-std'];
                        ColValues{end+1} = std(rdM(eve(:)==EVElabelID(j)));
                    end
                    
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-RD-mean'];
                        ColValues{end+1} = mean(rd1(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-RD-std'];
                        ColValues{end+1} = std(rd1(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-RD-mean'];
                        ColValues{end+1} = mean(rd2(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-RD-std'];
                        ColValues{end+1} = std(rd2(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-RD-mean'];
                        ColValues{end+1} = mean(rdM(bc(:)==BClabelID(j)));
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-RD-std'];
                        ColValues{end+1} = std(rdM(bc(:)==BClabelID(j)));
                    end
                    
                    %%%%%%%%%%%%%%% ROI Volumes
                    for j=1:length(BClabelNames)
                        ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'Volume'];
                        ColValues{end+1} = sum(bc(:)==BClabelID(j))*prod(bcinfo.PixelDimensions(1:3));
                    end
                    
                    for j=1:length(EVElabelNames)
                        ColHeader{end+1} = ['Eve-' EVElabelNames{j} '-' 'Volume'];
                        ColValues{end+1} = sum(eve(:)==EVElabelID(j))*prod(eveinfo.PixelDimensions(1:3));
                    end
                    
                    %%%%%%%%%%%%%%%  Write out
                    
                    if(~exist(AllStatsFile,'file'))
                        fp = fopen(AllStatsFile,'wt'); fprintf(fp,'%s, ',ColHeader{:}); fprintf(fp,'\n'); fclose(fp);
                    end
                    
                    fp = fopen(reportFileName,'at'); fprintf(fp,'%s, ',ColValues{1}); fprintf(fp,'%e, ',ColValues{2:end}); fprintf(fp,'\n'); fclose(fp);
                    fp = fopen(AllStatsFile,'at'); fprintf(fp,'%s, ',ColValues{1}); fprintf(fp,'%e, ',ColValues{2:end}); fprintf(fp,'\n'); fclose(fp);
                end
                %% Make a preview figure
                if(doMakePNG)
                    r = find(fa1(:,128,32));
                    x = round([min(r) median(r) max(r)]);
                    fa1(x,:,:)=1;
                    faM(x,:,:)=1;
                    eve(x,:,:)=128;
                    bc(x,:,:)=128;
                    figure(1);
                    clf;
                    imagesc(flipud([eve(:,:,32) faM(:,:,32)*255; bc(:,:,32) fa1(:,:,32)*255; fa2(:,:,32)*255 abs(fa1(:,:,32)-fa2(:,:,32))*255]') ); axis equal tight off
                    title([SUBJS(jSubj).name ' ' SESSIONS(jSession).name])
                    
                    drawnow
                    saveas(gcf,[D 'pngs' filesep SESSIONS(jSession).name '.png'])
                end
                
            end
        catch err
            fprintf('Oh no!\n')
            fp=fopen(FatalErrorFile,'at');
            fprintf(fp,'%s - %s (line %d)\n',SESSIONS(jSession).name,err.message, err(end).stack(end).line);
            fclose(fp);
        end
    end
end

cd(D2)