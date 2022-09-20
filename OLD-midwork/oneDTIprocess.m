% Setup the root directory for all subjects
%% ADDITIONAL SOFTWARE NEEDED
% fsl
% mrtrix

%% Original PATH
%D = '/fs4/masi/landmaba/BLSAdti/BLSA/'
%addpath('/fs4/masi/landmaba/BLSAdti/BLSA')
%addpath(genpath('/fs4/masi/landmaba/BLSAdti/matlab'))

%% Nazirah's Mac
% D = '/Users/nana/Documents/MATLAB/BLSA/'; % data path
% D2 = '/Users/nana/Documents/MATLAB/BLSA/parseXNAT_Nazirah/'; % output files, codes and labels path
% EVEpath = '/Users/nana/Documents/MATLAB/BLSA/';
% 
% addpath([D2 'Nazirah/functions/']);
% setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
% setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% addpath('/usr/local/bin/');

%% MASI Computers on Nazirah's Mac
% D = '/Users/nana/masi-42/Documents/data/';
% D2 = '/Users/nana/masi-42/Documents/Test1/';
% addpath([D2 'functions/']);

%% MASI Computers
%% MASI Computers - !!D and D2 must have filesep at the end!!
%D = '/home/local/VANDERBILT/mohdkhn/Documents/newdata/'; % data path
D = '/nfs2/harmonization/BLSA/'; %nfs2 datapath
D2 = '/home/local/VANDERBILT/mohdkhn/Documents/Test3-newdata/'; % code and labels path
EVEpath = '/nfs/masi/yangq6/EVE_Reg_BLSA/Reg/';

addpath([D2 'functions/']);
setenv('PATH', [getenv('PATH') ':/usr/local/anaconda3/bin']);
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

%% set time to record date and time to add to filename (FatalErrors.txt, AllStats.csv)
timedate = datestr(now,'yyyymmddTHHMMSS');
ReportFolder = [D2 'statsWithADRDVol' timedate];

if ~exist(ReportFolder,'dir')
    mkdir(ReportFolder)
end

FatalErrorFile = [D2 'FatalErrorsOneDTI-' timedate '.txt'];
AllStatsFile = [D2 'AllStats-HeaderWithADRDVol-' timedate '.csv'];

%%
% Load the label names
% EVE_path = [D 'EVE_Labels.csv'];
% BC_path = [D 'andrew_multiatlas_labels.csv'];
% DTI_path = [D 'DTIseg.txt'];

EVE1_path = [D2 'JHU_MNI_SS_WMPM_Type-I_SlicerLUT.txt'];
EVE2_path = [D2 'JHU_MNI_SS_WMPM_Type-II_SlicerLUT.txt'];
EVE3_path = [D2 'JHU_MNI_SS_WMPM_Type-III_SlicerLUT.txt'];
BC_path = [D2 'T1_label_volumes.txt'];
%DTI_path = [D2 'DTIseg.txt'];

%% Read labels
% EVE labels
EVE1labelNames = readcell(EVE1_path,"delimiter",'\n');
EVE1labelNames2 = split(EVE1labelNames(4:end,:),' ');
EVE1labelNames = join(EVE1labelNames2(:,[1 2]),' ');
for i=1:length(EVE1labelNames)
    EVE1labelID(i) = sscanf(EVE1labelNames{i},'%d');
end

EVE2labelNames = readcell(EVE2_path,"delimiter",'\n');
EVE2labelNames2 = split(EVE2labelNames(4:end,:),' ');
EVE2labelNames = join(EVE2labelNames2(:,[1 2]),' ');
for i=1:length(EVE2labelNames)
    EVE2labelID(i) = sscanf(EVE2labelNames{i},'%d');
end

EVE3labelNames = readcell(EVE3_path,"delimiter",'\n');
EVE3labelNames2 = split(EVE3labelNames(4:end,:),' ');
EVE3labelNames = join(EVE3labelNames2(:,[1 2]),' ');
for i=1:length(EVE3labelNames)
    EVE3labelID(i) = sscanf(EVE3labelNames{i},'%d');
end

% SLANT/Brain Color
BClabelNames = readcell(BC_path,"delimiter",'\n');
BClabelNames = BClabelNames(2:end); %remove header
BClabelNames2 = split(BClabelNames,',');
BClabelNames = join(BClabelNames2(:,[2 1]),' ');
for i=1:length(BClabelNames2)
    BClabelID(i) = str2double(BClabelNames2{i,2});
end

%[EVElabelNames,EVElabelID,BClabelNames,BClabelID,DTISeglabelNames,DTISeglabelID] = get_label_names(EVE_path,BC_path,DTI_path);

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

for jSession=1:length(SESSIONS)
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
    test = 0;
    try
        % RULE:
        % dtiQA must be more 1
        % WM stamper must exist (equals to 1)
        % Slant must exist (equals to 1)
        
        % IF dtiQA < 2 and HAS Stamper AND MultiAtlas
        if test ==1 %~(length(dtiQA)>1 && length(Stamper)==1 && length(Slant)==1) && length(dtiQA) == 1
            % The data are not there. let's write a stats file that
            % tells the user why.
            
            %if length(dtiQA)<2
                % This write report only if dtiQA is less than 2.
                % Stamper and MultiAtlas must exist.
                %reportFileName = [ReportFolder filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                %fp = fopen(reportFileName,'at');
                %fprintf(fp,'ONLY 1 DTI');
                %fclose(fp);
                %error('There are not 2 DTIs.');
            %end
            
            
        %NKcomment else % Now, we have enough DTI files, let's check Stamper, CSV, and PNG files
            % if(and(and(length(dtiQA)>1,length(Stamper)==1),length(MultiAtlas)==1))
            
            % CHECK: if WM Stamper folder has less than 3 files, don't continue. Move to the next session.
            %if(length(dir([DS Stamper(1).name]))<3)
            %    error('WM Stamper has less than 3 files')
            %end
            fprintf('FOUND valid dataset\n')
            
            % Create CSV report file.
            %reportFileName = [ReportFolder filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
            
            % If CSV file exist, don't continue below. Move to the next Session
            % Nazirah: this checker is useful if change reportFolder to no-date name
            %if(exist(reportFileName,'file'))
            %    fprintf(['Already done: ' reportFileName '\n']);
            %    continue;
            %end
            %fprintf('touch\n')
%             system(['echo `date` > ' reportFileName]);
            
            % If PNG exists, don't continue below and go to next session. [Whyy?]
%             if(doSkipPNGDone)
%                 if(exist([D 'pngs' filesep SESSIONS(jSession).name '.png'],'file'))
%                     fprintf('It''s done\n');
%                     continue;
%                 end
%             end
            
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
            
            %% Deal with the first DTI session "DTI(1)"
            %dti1
            DTI1 = 1;
            cd([DS dtiQA(DTI1).name filesep 'SCALARS'])
            fa1name = findfileniiorgz([pwd filesep], 'fa.nii');
            md1name= findfileniiorgz([pwd filesep], 'md.nii');
            ad1name = findfileniiorgz([pwd filesep], 'ad.nii');
            rd1name = findfileniiorgz([pwd filesep], 'rd.nii' );
            
            
            %% Resample the EVE and BrainColor masks from T1 to DTI(1) space
            
            % Use flirt to find a 6 dof transform
            %                 xfmname = [DS 'mpr2fa.txt'];
            %                 system(['flirt -in ' mprname ' -out ' mprname '-flirt.nii.gz -ref ' faMname ' -omat ' xfmname ' -dof 6'])
            %                 system(['flirt -in ' label1name ' -ref ' faMname ' -applyxfm -init ' xfmname ' -interp nearestneighbour' ' -out ' label1name '-flirt.nii.gz'])
            %                 system(['flirt -in ' label1name ' -ref ' faMname ' -applyxfm -init ' xfmname ' -interp nearestneighbour'])
            
            % check if xfm matrix is there: WE NOW USE "mprage2b0.mat"
            xfm1name = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'mprage2b0.mat'];
            
            % if xfm not available, need to find the xfm first!
            if(length(dir(xfm1name))<1)
                
                mprbrainname = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'mpr_brain.nii.gz'];
                b0name = [DS dtiQA(DTI1).name filesep 'PREPROCESSED' filesep 'b0.nii.gz'];
                dwi1name = [DS dtiQA(DTI1).name filesep 'PREPROCESSED' filesep 'dwmri.nii.gz'];
                
                
                % if b0 is not there, get the b0 first
                if (length(dir(b0name))<1)
                    bvec1name = [DS dtiQA(DTI1).name filesep 'PREPROCESSED' filesep 'dwmri.bvec'];
                    bval1name = [DS dtiQA(DTI1).name filesep 'PREPROCESSED' filesep 'dwmri.bval'];
                    system(['dwiextract ' dwi1name ' - -bzero -fslgrad ' bvec1name ' ' bval1name ' -quiet | mrmath - mean ' b0name ' -axis 3 -quiet'])
                end
                    
                
                % if MPRAGE brain is not there, need to get the skull stripped first
                if length(dir(mprbrainname))<1
                    system(['bet ' mprname ' ' mprbrainname ' -R -f 0.5 -g 0 -m']);
                end
                
                % now we have b0, mprbrain --> find the xfm1
                xfm1invname = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'b02mprage'];
                system(['epi_reg --epi=' b0name ' --t1=' mprname ' --t1brain=' mprbrainname ' --out=' xfm1invname]);
                system(['convert_xfm -omat ' xfm1name ' -inverse ' xfm1invname '.mat']);
                
                % 
                %dwi1file = dir([DS2 filesep 'DTI1' filesep 'NIFTI' filesep '*DTI1.nii.gz']);
                %dwi1name = [dwi1file(1).folder filesep dwi1file(1).name];
                %dwi2mpragename = [dwi1file(1).folder filesep 'dwi2mprage.nii.gz'];
                %system(['flirt -in ' mprbrainname ' -out ' mprbrainname '-flirt.nii.gz -ref ' fa1name ' -omat ' xfm1name ' -dof 6']); 
                %system(['epi_reg --epi=' dwi1name ' --t1=' mprname ' --t1brain=' mprbrainname ' --out=' dwi2mpragename]) 
            end
            
            % Resample the FA labels from Eve into Subject FA
%             xfm1name = [DS Stamper(1).name filesep '/TRANSFORMATION/fa_t1_transformation_matrix.txt'];
%             if(length(dir(xfm1name))<1)
%                 xfm1name = [DS Stamper(1).name filesep '/Intra_Session_Reg/outputAffine.txt'];
%             end
%             if(length(dir(xfm1name))<1)
%                 error('Cannot find WM transform');
%             end
            label1name = [DS Stamper(1).name filesep 'WM_LABELS' filesep 'EVE1_Labels.nii.gz'];
            label2name = [DS Stamper(1).name filesep 'WM_LABELS' filesep 'EVE2_Labels.nii.gz'];
            label3name = [DS Stamper(1).name filesep 'WM_LABELS' filesep 'EVE3_Labels.nii.gz'];
    
            eve1Name = [label1name '.subjLabels.nii.gz'];
            eve2Name = [label2name '.subjLabels.nii.gz'];
            eve3Name = [label3name '.subjLabels.nii.gz'];
            
            if length(dir(eve1Name))<1 || length(dir(eve2Name))<1 || length(dir(eve3Name))<1
                
                if length(dir(label1name))<1 || length(dir(label2name))<1 || length(dir(label3name))<1
                    % if labels doesn't exist in WM_LABELS, copy from Qi's EVEpath
                    copyfile([EVEpath SESSIONS(jSession).name filesep '*Type-I.nii.gz'],label1name)
                    copyfile([EVEpath SESSIONS(jSession).name filesep '*Type-II.nii.gz'],label1name)
                    copyfile([EVEpath SESSIONS(jSession).name filesep '*Type-III.nii.gz'],label1name)
                end
                    
                % resample EVE to FA (use dwi2mprage to inverse to mprage2dwi then register)
                %system(['convert_xfm -omat ' xfm1name '.inv -inverse' xfm1name])
                system(['flirt -in ' label1name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eve1Name]);
                system(['flirt -in ' label2name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eve2Name]);
                system(['flirt -in ' label3name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eve3Name]);
            end
            
            
            % Resample the Multi-Atlas Labels
            masegfile = dir([DS Slant(1).name filesep 'SEG' filesep 'T1_seg.nii*']);
            masegname = [DS Slant(1).name filesep 'SEG' filesep masegfile(1).name];
            
            brainColorName = [masegname '.subjLabels.nii.gz'];
            if(length(dir(brainColorName))<1)
                system(['flirt -in ' masegname ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' brainColorName]); 
                %system(['reg_resample -aff ' xfm1name '.inv ' '-ref ' fa1name ' -flo ' masegname ' -res ' masegname '.subjLabels.nii.gz' ' -inter 0'])
            end
        elseif ~(length(dtiQA)>1 && length(dtiQAMulti)==1 && length(Stamper)==1 && length(Slant)==1)
            error('dtiQA = %d, dtiQA_double = %d, WMStamper = %d, Slant = %d, MPRAGE = %d',length(dtiQA), length(dtiQAMulti), length(Stamper), length(Slant), length(MPRAGE))
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