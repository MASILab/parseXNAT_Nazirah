% Setup the root directory for all subjects

%% Original PATH
%D = '/fs4/masi/landmaba/BLSAdti/BLSA/'
%addpath('/fs4/masi/landmaba/BLSAdti/BLSA')
%addpath(genpath('/fs4/masi/landmaba/BLSAdti/matlab'))

%% Nazirah's Mac
D = '/Users/nana/Documents/MATLAB/BLSA/';
D2 = '/Users/nana/Documents/MATLAB/BLSA/parseXNAT_Nazirah/';

addpath([D2 'Nazirah/functions/']);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
%addpath('/usr/local/bin/');

%% MASI Computers on Nazirah's Mac
% D = '/Users/nana/masi-42/Documents/data/';
% D2 = '/Users/nana/masi-42/Documents/Test1/';
% addpath([D2 'functions/']);

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
BC_path = [D2 'T1_label_volumes.txt'];
DTI_path = [D2 'DTIseg.txt'];

[EVElabelNames,EVElabelID,BClabelNames,BClabelID,DTISeglabelNames,DTISeglabelID] = get_label_names(EVE_path,BC_path,DTI_path);

%% Get a list of Subjects; ignore . and ..
%SUBJS = dir([D filesep 'BLSA*']); SUBJS=SUBJS(1:end);



%% This is used to ensure that the stats file is replaced each time.
doMakePNG = 0;
doSkipPNGDone = 0;
doMakeReport = 1;


%% Let's get the stats file ready!

% Get a list of Sessions; ignore . and ..
SESSIONS = dir([D 'BLSA*']);
SESSIONS(startsWith({SESSIONS.name},'.')) = [];

for jSession=2%1:length(SESSIONS)
    fprintf('%s\n', SESSIONS(jSession).name)
    
    DS = [D SESSIONS(jSession).name filesep 'ASSESSORS' filesep];
    fprintf('DS: %s\n',DS);
    DS2 = [D SESSIONS(jSession).name filesep 'SCANS' filesep];
    fprintf('DS2: %s\n',DS2);
    
    % Find Stamper folder
    Stamper = dir([DS  '*Stamper*']);
    
    % Find all dtiQA folder - focus on v7
    dtiQA = dir([DS '*dtiQA_synb0_v7*']);
    if (length(dtiQA))>1
        fprintf('Got %d dtiQA_v7\n',length(dtiQA))
    end
    
    % Find dtiQAMulti - dtiQA*double*
    dtiQAMulti = dir([DS '*dtiQA*double*']);
    
    % Find Slant folder (replace MPRAGE*Multi_Atlas)
    Slant = dir([DS '*slant*']);
    
    % Sometimes Multi-atlas will find another T1. If
    % more than one are found, choose the MPRAGE one.
    
    % MPRAGE is located in SCANS folder
    MPRAGE = dir([DS2 'MPRAGE*']);
    
    % Checking if all is there
    fprintf([SESSIONS(jSession).name '\n']);
    fprintf('Stamper: %d dtiQA: %d Slant: %d dtiQAMulti: %d\n', length(Stamper),length(dtiQA), length(Slant), length(dtiQAMulti));
    
    try
        % RULE:
        % dtiQA must be more 1
        % WM stamper must exist (equals to 1)
        % Slant must exist (equals to 1)
        
        % IF dtiQA < 2 and HAS Stamper AND MultiAtlas
        if ~(length(dtiQA)>1 && length(Stamper)==1 && length(Slant)==1)
            % The data are not there. let's write a stats file that
            % tells the user why.
            
            if length(dtiQA)<2
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
            % Location: SESSION_NAME > SCANS > NIFTI > filename.gz
            mprfile = dir([DS2 MPRAGE(1).name filesep 'NIFTI' filesep '*MPRAGE*.gz']);
            mprname = [mprfile(1).folder filesep mprfile(1).name];
            %mprname = [DS MPRAGE(1).name filesep 'NIFTI' filesep mprfile(1).name];
            
            
            %% Deal with the Multi DTI session "DTI multi"/"DTI double"
            % Location: SESSION_NAME > ASSESSORS > dtiQAMulti > SCALARS
            cd([DS dtiQAMulti(1).name filesep 'SCALARS'])
            adMname = findfileniiorgz([pwd filesep], 'ad.nii');
            rdMname = findfileniiorgz([pwd filesep], 'rd.nii');
            faMname = findfileniiorgz([pwd filesep],'fa.nii');
            mdMname = findfileniiorgz([pwd filesep],'md.nii');
            %roiMname = findfileniiorgz([pwd filesep 'extra'], 'multi_atlas_labels.nii');
            
            %% find which is DT1 and DTI2
            %DTI1 = 1;
            %DTI2 = 2;
            for i = 1:length(dtiQA)
                outlogfile = dir([DS dtiQA(i).name filesep 'OUTLOG' filesep '*.txt']);
                outlog = fileread([outlogfile.folder filesep outlogfile.name]);
                
                foundDTI1 = strfind(outlog,'DTI1');
                foundDTI2 = strfind(outlog,'DTI2');
                
                if ~isempty(foundDTI1) && isempty(foundDTI2) % found DT1
                    DTI1 = i;
                elseif isempty(foundDTI1) && ~isempty(foundDTI2) % found DT2
                    DTI2 = i;
                else
                    error('DTI1 and DTI2 not found or both in the same directory')
                end
            end
            
            
            %% Deal with the first DTI session "DTI(1)"
            %dti1
            cd([DS dtiQA(DTI1).name filesep 'SCALARS'])
            fa1name = findfileniiorgz([pwd filesep], 'fa.nii');
            md1name= findfileniiorgz([pwd filesep], 'md.nii');
            ad1name = findfileniiorgz([pwd filesep], 'ad.nii');
            rd1name = findfileniiorgz([pwd filesep], 'rd.nii' );
            
            roi1name = [pwd filesep 'multi_atlas_labels.nii.gz'];
            
            %% Resample the EVE and BrainColor masks from T1 to DTI(1) space
            masegfile = dir([DS Slant(1).name filesep 'SEG' filesep 'T1_seg.nii*']);
            masegname = [DS Slant(1).name filesep 'SEG' filesep masegfile(1).name];
            
            % Use flirt to find a 6 dof transform
            %                 xfmname = [DS 'mpr2fa.txt'];
            %                 system(['flirt -in ' mprname ' -out ' mprname '-flirt.nii.gz -ref ' faMname ' -omat ' xfmname ' -dof 6'])
            %                 system(['flirt -in ' label1name ' -ref ' faMname ' -applyxfm -init ' xfmname ' -interp nearestneighbour' ' -out ' label1name '-flirt.nii.gz'])
            %                 system(['flirt -in ' label1name ' -ref ' faMname ' -applyxfm -init ' xfmname ' -interp nearestneighbour'])
            
            % check if xfm matrix is there: WE NOW USE "mpr2fa_transformation_matrix.txt"
            xfm1name = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'mpr2fa.txt'];
            if(length(dir(xfm1name))<1)
                % if not available, run flirt to find xfm matrix
                mprbrainname = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'mpr_brain.nii.gz'];
                
                if length(dir(mprbrainname))<1
                    system(['bet ' mprname ' ' mprbrainname ' -R -f 0.5 -g 0 -m']);
                end
                
                system(['flirt -in ' mprbrainname ' -out ' mprbrainname '-flirt.nii.gz -ref ' fa1name ' -omat ' xfm1name ' -dof 6']);             
            end
            
            % Resample the FA labels from Eve into Subject FA
%             xfm1name = [DS Stamper(1).name filesep '/TRANSFORMATION/fa_t1_transformation_matrix.txt'];
%             if(length(dir(xfm1name))<1)
%                 xfm1name = [DS Stamper(1).name filesep '/Intra_Session_Reg/outputAffine.txt'];
%             end
%             if(length(dir(xfm1name))<1)
%                 error('Cannot find WM transform');
%             end
            label1name = [DS Stamper(1).name filesep 'WM_LABELS' '/Rectified_EVE_Labels.nii.gz'];
            eveName = [label1name '.subjLabels.nii.gz'];
            %if(length(dir(eveName))<1)
                
                % resample EVE to FA
                system(['flirt -in ' label1name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eveName]);
                
                %% REMOVED
                % reg_transform -ref ref_img -invAff transform_mat transform_mat.inv
                % link: http://cmictig.cs.ucl.ac.uk/wiki/index.php/Reg_transform
                % system(['reg_transform -ref ' fa1name ' -invAffine ' xfm1name ' ' xfm1name '.inv']);
                % input: xfm1name || output: xfm1name.inv
                
                % reg_resample -aff transform_mat.inv -ref ref_img -flo EVE_Labels.nii.gz -res EVE_Labels.nii.gz.subjLabels.nii.gz -inter 0
                % link: http://cmictig.cs.ucl.ac.uk/wiki/index.php/Reg_resample
                % input: xfm1name, label1name, fa1name || output: label1name.subjLabels.nii.gz
                % system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' label1name ' -res ' label1name '.subjLabels.nii.gz' ' -inter 0'])
                %system(['flirt -in ' label1name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' label1name '.subjLabels.nii.gz'])
            %end
            
            
            % Resample the Multi-Atlas Labels
            brainColorName = [masegname '.subjLabels.nii.gz'];
            %if(length(dir(brainColorName))<1)
                system(['flirt -in ' masegname ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' brainColorName]); 
                %system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' masegname ' -res ' masegname '.subjLabels.nii.gz' ' -inter 0'])
            %end
            
            
            %% Now deal with the second DTI session
            
            cd([DS dtiQA(DTI2).name filesep 'SCALARS'])
            ad2name = findfileniiorgz([pwd filesep], 'ad.nii');
            rd2name = findfileniiorgz([pwd filesep], 'rd.nii');
            fa2name = findfileniiorgz([pwd filesep],'fa.nii');
            md2name = findfileniiorgz([pwd filesep],'md.nii');
            %roi2name = findfileniiorgz([pwd filesep], 'multi_atlas_labels.nii');
            
            %% Verify
            
            %faM = loaduntouchniiorniigz(faMname);
            %faM2 = niftiread(faMname);
            
            faM = loadniiorgz(faMname); faMinfo = infoniiorgz(faMname);
            fa1 = loadniiorgz(fa1name); fa1info = infoniiorgz(fa1name);
            fa2 = loadniiorgz(fa2name); fa2info = infoniiorgz(fa2name);
            
            mdM = loadniiorgz(mdMname);
            md1 = loadniiorgz(md1name);
            md2 = loadniiorgz(md2name);
            
            adM = loadniiorgz(adMname);
            ad1 = loadniiorgz(ad1name);
            ad2 = loadniiorgz(ad2name);
            
            rdM = loadniiorgz(rdMname);
            rd1 = loadniiorgz(rd1name);
            rd2 = loadniiorgz(rd2name);
            
            eve = loadniiorgz(eveName); eveinfo = infoniiorgz(eveName);
            bc = loadniiorgz(brainColorName); bcinfo = infoniiorgz(brainColorName);
            
            %roi1 = loadniiorgz(roi1name);
            %roi2 = loadniiorgz(roi2name);
            
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
                
                %% PG 6 - DSC on labels
%                 for j=1:length(DTISeglabelNames)
%                     ColHeader{end+1} = ['STATS-ROI-' DTISeglabelNames{j} '-DSC'];
%                     A = roi1==j;
%                     B = roi2==j;
%                     ColValues{end+1} = 2 * sum(A(:).*B(:)) / (sum(A(:))+sum(B(:)));
%                 end
                
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
        fprintf('OH NOOOOO!\n')
        fprintf('%s - %s (line %d)\n',SESSIONS(jSession).name,err.message, err(end).stack(end).line);
        fp=fopen(FatalErrorFile,'at');
        fprintf(fp,'%s - %s (line %d)\n',SESSIONS(jSession).name,err.message, err(end).stack(end).line);
        fclose(fp);
    end
end

cd(D2)