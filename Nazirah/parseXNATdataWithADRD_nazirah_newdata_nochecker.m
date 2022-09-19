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
D2 = '/home/local/VANDERBILT/mohdkhn/Documents/Test4-newdata/'; % code and labels path
EVEpath = '/nfs/masi/yangq6/EVE_Reg_BLSA/*/';

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

FatalErrorFile = [D2 'FatalErrorsWithADRDv8-' timedate '.txt'];
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
% SESSIONS = dir([D 'BLSA_*']);
% SESSIONS(startsWith({SESSIONS.name},'.')) = [];
SESSIONS_LIST = readcell('/home/local/VANDERBILT/mohdkhn/Documents/Test4-newdata/IncompleteDirectories.xlsx');
SESSIONS_LIST = SESSIONS_LIST(strcmp(SESSIONS_LIST(:,17),'Include'),1);
SESSIONS = cell2struct(SESSIONS_LIST(:,1)','name',1);
fprintf('Total sessions to run: %d \n',length(SESSIONS));

for jSession=1:length(SESSIONS)
    
    rerun = 0;
    fprintf('\n\n=============\n')
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
    dtiQAMulti = dir([DS '*dtiQA*double_v7*']);
    
    % Find Slant folder - *slant*
    Slant = dir([DS '*slant*']);
    
    % MPRAGE is located in SCANS folder
    MPRAGE = dir([DS2 'MPRAGE*']);
    
    % Checking if all is there
    %fprintf([SESSIONS(jSession).name '\n']);
    fprintf('--> Stamper=%d dtiQA=%d Slant=%d dtiQAMulti=%d MPRAGE=%d \n', length(Stamper),length(dtiQA), length(Slant), length(dtiQAMulti), length(MPRAGE));
    
    try
        % RULE:
        % dtiQA_synb0_v7 must be >= 1. If 1, one stats must be NaNs.
        % dtiQA_synb0_double_v7 must be = 1, if not, check which one
        % WM stamper must exist (equals to 1)
        % Slant must exist (equals to 1). If >=1, check which one.
        % MPRAGE must exist (equals to 1)
        
        % CHECK1: check if dirs not exist at all --> send to ERROR
        dircounts = [length(Stamper) length(dtiQA) length(Slant) length(dtiQAMulti) length(MPRAGE)];
        if dircounts(4)==0 %only dtiQAMulti=0
            err_message = 'No dtiQAMulti. DTIM stats will be NaNs';
            write_warning_message(FatalErrorFile, SESSIONS(jSession).name, err_message);
        elseif sum(dircounts==0)>0
            error('Stamper=%d dtiQA=%d Slant=%d dtiQAMulti=%d MPRAGE=%d ', dircounts(1),dircounts(2),dircounts(3),dircounts(4),dircounts(5));
        end
        
        % CHECK2: check Slant --> choose which Slant
        if length(Slant)>1
            % if more than 2, choose slant from MPRAGE, not other T1
            for i = 1:length(Slant)
                outlogfile = dir([DS Slant(i).name filesep 'PBS' filesep '*.slurm']);
                outlog = fileread([outlogfile.folder filesep outlogfile.name]);
                
                foundSlantMPRAGE = strfind(outlog,'MPRAGE');
                
                if foundSlantMPRAGE
                    slantID = i;
                end
            end
            
            if ~exist('slantID','var')
                error('Slant exist but could not find correct Slant.')
            end
            
            %overwrite Slant to only at slantID
            Slant = Slant(slantID);
            
            % put a warning message in FatalErrors
            err_message = sprintf('Slant used is %s',Slant.name);
            write_warning_message(FatalErrorFile, SESSIONS(jSession).name, err_message);
        end
        
        % CHECK3: check dtiQA_double --> choose which dtiQA_double
        % Main pref: dtiQA_synb0_double
        % if no synb0 version, choose the regular dtiQA_double
        % if there is 1 synb0
        if length(dtiQAMulti)>1
            % if more than 2, choose slant from MPRAGE, not other T1
            
            for i = 1:length(dtiQAMulti)
                
                if contains(dtiQAMulti(i).name, 'synb0')
                    dtiQAMultiID = i;
                    break;
                else
                    outlogfile = dir([DS dtiQAMulti(i).name filesep 'PBS' filesep '*.slurm']);
                    outlog = fileread([outlogfile.folder filesep outlogfile.name]);
                    
                    foundT1WMT = strfind(outlog,'T1WMT');
                    foundT1SPGR = strfind(outlog,'T1_SPGR');
                    
                    if isempty(foundT1WMT) && isempty(foundT1SPGR)
                        dtiQAMultiID = i;
                    end
                end
            end
            
            if ~exist('dtiQAMultiID','var')
                error('dtiQA_synb0_double found but no with only [DTI1,DTI2,MPRAGE]')
            end
            
            % overwrite dtiQAMulti to only at dtiQAMultiID
            dtiQAMulti = dtiQAMulti(dtiQAMultiID);
            
            % put a warning message in FatalErrors
            err_message = sprintf('dtiQA_double used is %s',dtiQAMulti.name);
            write_warning_message(FatalErrorFile, SESSIONS(jSession).name, err_message);
        end
        
        % Note:
        % At this line, dtiQAMulti and Slant should equals 1 only.
        % Warning: dtiQAMulti = 0 --> will deal with it below
        % Now will check how many dtiQA (should have 1 or 2)
        
        % IF have all directories
        if ~(~isempty(dtiQA) && length(Stamper)==1 && length(Slant)==1)
            
            if length(dtiQA)<1
                % This write report only if dtiQA is less 1
                reportFileName = [ReportFolder filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                fp = fopen(reportFileName,'at');
                fprintf(fp,'NO DTI');
                fclose(fp);
            end
            
            
        else % Now, we have enough DTI files, let's check Stamper, CSV, and PNG files
            % if(and(and(length(dtiQA)>1,length(Stamper)==1),length(MultiAtlas)==1))
            
            % CHECK: if WM Stamper folder has less than 3 files, don't continue. Move to the next session.
            if(length(dir([DS Stamper(1).name]))<3)
                %continue;
            end
            fprintf('FOUND valid dataset\n')
            
            % Create CSV report file.
            reportFileName = [ReportFolder filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
            
            % If CSV file exist, don't continue below. Move to the next Session
            % Nazirah: this checker is useful if change reportFolder to no-date name
            if(exist(reportFileName,'file'))
                fprintf(['Already done: ' reportFileName '\n']);
                continue;
            end
            
            %% NOW, ALL FILES EXISTS AND NEED TO BE RUN
            %% Find the MPRAGE
            % Location: SESSION_NAME > SCANS > NIFTI > filename.gz
            mprfile = dir([DS2 MPRAGE(1).name filesep 'NIFTI' filesep '*MPRAGE*.gz']);
            mprname = [mprfile(1).folder filesep mprfile(1).name];
            %mprname = [DS MPRAGE(1).name filesep 'NIFTI' filesep mprfile(1).name];
            
            
            %% Deal with the Multi DTI session "DTI multi"/"DTI double"
            % Location: SESSION_NAME > ASSESSORS > dtiQAMulti > SCALARS
            if ~isempty(dtiQAMulti)
                DTIM = 1;
                cd([DS dtiQAMulti(1).name filesep 'SCALARS'])
                adMname = findfileniiorgz([pwd filesep], 'ad.nii');
                rdMname = findfileniiorgz([pwd filesep], 'rd.nii');
                faMname = findfileniiorgz([pwd filesep],'fa.nii');
                mdMname = findfileniiorgz([pwd filesep],'md.nii');
            elseif isempty(dtiQAMulti)
                DTIM = 0;
            end
            
            %roiMname = findfileniiorgz([pwd filesep 'extra'], 'multi_atlas_labels.nii');
            
            %% find which is DT1 and DTI2
            DTI1 = 0;
            DTI2 = 0;
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
            
            if DTI1 == 0 || DTI2 == 0
                err_message = 'There is only 1 DTI. The other DTI stats will be NaNs';
                write_warning_message(FatalErrorFile, SESSIONS(jSession).name, err_message);
            end
            
            %% Deal with the first DTI session "DTI(1)"
            if DTI1 > 0
                %dti1
                cd([DS dtiQA(DTI1).name filesep 'SCALARS'])
                fa1name = findfileniiorgz([pwd filesep], 'fa.nii');
                md1name= findfileniiorgz([pwd filesep], 'md.nii');
                ad1name = findfileniiorgz([pwd filesep], 'ad.nii');
                rd1name = findfileniiorgz([pwd filesep], 'rd.nii' );
            end
            
            %% Now deal with the second DTI session
            if DTI2 > 0
                cd([DS dtiQA(DTI2).name filesep 'SCALARS'])
                ad2name = findfileniiorgz([pwd filesep], 'ad.nii');
                rd2name = findfileniiorgz([pwd filesep], 'rd.nii');
                fa2name = findfileniiorgz([pwd filesep],'fa.nii');
                md2name = findfileniiorgz([pwd filesep],'md.nii');
            end
            
            %% Resample the EVE and BrainColor masks from T1 to DTI(1) space
            if DTI1 > 0
                dtireg = DTI1;
            elseif DTI2 > 0
                dtireg = DTI2;
            end
            
            % check if xfm matrix is there: WE NOW USE "mprage2b0.mat"
            xfm1name = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'mprage2b0.mat'];
            
            % if xfm not available, need to find the xfm first!
            if(length(dir(xfm1name))<1)
                
                mprbrainname = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'mpr_brain.nii.gz'];
                b0name = [DS dtiQA(dtireg).name filesep 'PREPROCESSED' filesep 'b0.nii.gz'];
                dwi1name = [DS dtiQA(dtireg).name filesep 'PREPROCESSED' filesep 'dwmri.nii.gz'];
                
                
                % if b0 is not there, get the b0 first
                if (length(dir(b0name))<1)
                    bvec1name = [DS dtiQA(dtireg).name filesep 'PREPROCESSED' filesep 'dwmri.bvec'];
                    bval1name = [DS dtiQA(dtireg).name filesep 'PREPROCESSED' filesep 'dwmri.bval'];
                    system(['dwiextract ' dwi1name ' - -bzero -fslgrad ' bvec1name ' ' bval1name ' -quiet | mrmath - mean ' b0name ' -axis 3 -quiet']);
                end
                
                
                % if MPRAGE brain is not there, need to get the skull stripped first
                if length(dir(mprbrainname))<1
                    system(['bet ' mprname ' ' mprbrainname ' -R -f 0.5 -g 0 -m']);
                end
                
                % now we have b0, mprbrain --> find the xfm1
                xfm1invname = [DS2 MPRAGE(1).name filesep 'NIFTI' filesep 'b02mprage'];
                system(['epi_reg --epi=' b0name ' --t1=' mprname ' --t1brain=' mprbrainname ' --out=' xfm1invname]);
                system(['convert_xfm -omat ' xfm1name ' -inverse ' xfm1invname '.mat']); % convert b02mprage to mprage2b0
                
            end
            
            % Resample the FA labels from Eve into Subject FA
            label1name = [DS Stamper(1).name filesep 'WM_LABELS' filesep 'EVE1_Labels.nii.gz'];
            label2name = [DS Stamper(1).name filesep 'WM_LABELS' filesep 'EVE2_Labels.nii.gz'];
            label3name = [DS Stamper(1).name filesep 'WM_LABELS' filesep 'EVE3_Labels.nii.gz'];
            
            eve1Name = [label1name '.subjLabels.nii.gz'];
            eve2Name = [label2name '.subjLabels.nii.gz'];
            eve3Name = [label3name '.subjLabels.nii.gz'];
            
            if length(dir(eve1Name))<1 || length(dir(eve2Name))<1 || length(dir(eve3Name))<1
                % if one of the labels are not there, need to register (create .subjLabels.nii.gz files)
                
                if length(dir(label1name))<1 || length(dir(label2name))<1 || length(dir(label3name))<1
                    % if labels doesn't exist in WM_LABELS, copy from Qi's EVEpath
                    system(['cp ' EVEpath SESSIONS(jSession).name filesep '*Type-I.nii.gz ' label1name]);
                    system(['cp ' EVEpath SESSIONS(jSession).name filesep '*Type-II.nii.gz ' label2name]);
                    system(['cp ' EVEpath SESSIONS(jSession).name filesep '*Type-III.nii.gz ' label3name]);
                    %copyfile([EVEpath SESSIONS(jSession).name filesep '*Type-I.nii.gz'],label1name)
                    %copyfile([EVEpath SESSIONS(jSession).name filesep '*Type-II.nii.gz'],label2name)
                    %copyfile([EVEpath SESSIONS(jSession).name filesep '*Type-III.nii.gz'],label3name)
                end
                
                % resample all EVEs to FA using mprage2b0
                system(['flirt -in ' label1name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eve1Name]);
                system(['flirt -in ' label2name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eve2Name]);
                system(['flirt -in ' label3name ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' eve3Name]);
            end
            
            
            % Resample the Multi-Atlas Labels
            masegfile = dir([DS Slant(1).name filesep 'SEG' filesep 'T1_seg.nii*']);
            masegname = [DS Slant(1).name filesep 'SEG' filesep masegfile(1).name];
            
            brainColorName = [masegname '.subjLabels.nii.gz'];
            if(length(dir(brainColorName))<1)
                % if subjLabels.nii.gz not there, need to register
                system(['flirt -in ' masegname ' -ref ' fa1name ' -applyxfm -init ' xfm1name ' -interp nearestneighbour' ' -out ' brainColorName]);
            end
            
            
            %% Verify
            
            % dti1
            if DTI1 > 0
                fa1 = loadniiorgz(fa1name); fa1info = infoniiorgz(fa1name);
                md1 = loadniiorgz(md1name);
                ad1 = loadniiorgz(ad1name);
                rd1 = loadniiorgz(rd1name);
            end
            
            % dti2
            if DTI2 > 0
                fa2 = loadniiorgz(fa2name); fa2info = infoniiorgz(fa2name);
                md2 = loadniiorgz(md2name);
                ad2 = loadniiorgz(ad2name);
                rd2 = loadniiorgz(rd2name);
            end
            
            % dtiM
            if DTIM > 0
                faM = loadniiorgz(faMname); faMinfo = infoniiorgz(faMname);
                mdM = loadniiorgz(mdMname);
                adM = loadniiorgz(adMname);
                rdM = loadniiorgz(rdMname);
            end
            
            % atlases
            eve1 = loadniiorgz(eve1Name); eve1info = infoniiorgz(eve1Name);
            eve2 = loadniiorgz(eve2Name); eve2info = infoniiorgz(eve2Name);
            eve3 = loadniiorgz(eve3Name); eve3info = infoniiorgz(eve3Name);
            bc = loadniiorgz(brainColorName); bcinfo = infoniiorgz(brainColorName);
            
            % save headercheck
            fp = fopen([D2 'headerCheckv8.txt'],'at');
            fprintf(fp,'%s ',SESSIONS(jSession).name);
            fprintf(fp,'%f %f %f %f %f %f %f %f %f %f %f %f ', ...
                bcinfo.PixelDimensions(1),bcinfo.PixelDimensions(2),bcinfo.PixelDimensions(3),...
                eve1info.PixelDimensions(1),eve1info.PixelDimensions(2),eve1info.PixelDimensions(3),...
                eve2info.PixelDimensions(1),eve2info.PixelDimensions(2),eve2info.PixelDimensions(3),...
                eve3info.PixelDimensions(1),eve3info.PixelDimensions(2),eve3info.PixelDimensions(3));
            
            if DTI1 > 0
                fprintf(fp,'%f %f %f ', fa1info.PixelDimensions(1),fa1info.PixelDimensions(2),fa1info.PixelDimensions(3));
            elseif DTI1 == 0
                fprintf(fp,'%f %f %f ',NaN,NaN,NaN);
            end
            
            if DTI2 > 0
                fprintf(fp,'%f %f %f ', fa2info.PixelDimensions(1),fa2info.PixelDimensions(2),fa2info.PixelDimensions(3));
            elseif DTI2 == 0
                fprintf(fp,'%f %f %f ',NaN,NaN,NaN);
            end
            
            if DTIM > 0
                fprintf(fp,'%f %f %f ', faMinfo.PixelDimensions(1),faMinfo.PixelDimensions(2),faMinfo.PixelDimensions(3));
            elseif DTIM == 0
                fprintf(fp,'%f %f %f ',NaN,NaN,NaN);
            end
            
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
                % EVE1
                for j=1:length(EVE1labelNames)
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-FA-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-FA-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(fa1(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(fa1(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-FA-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-FA-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(fa2(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(fa2(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-FA-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-FA-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(faM(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(faM(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                % EVE2
                for j=1:length(EVE2labelNames)
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-FA-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-FA-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(fa1(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(fa1(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-FA-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-FA-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(fa2(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(fa2(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-FA-mean'];
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-FA-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(faM(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(faM(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                % EVE3
                for j=1:length(EVE3labelNames)
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-FA-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-FA-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(fa1(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(fa1(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-FA-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-FA-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(fa2(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(fa2(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-FA-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-FA-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(faM(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(faM(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                % SLANT BC
                for j=1:length(BClabelNames)
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-FA-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-FA-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(fa1(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(fa1(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-FA-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-FA-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(fa2(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(fa2(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-FA-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-FA-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(faM(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(faM(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                % MD
                % EVE1
                for j=1:length(EVE1labelNames)
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-MD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-MD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(md1(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(md1(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-MD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-MD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(md2(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(md2(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-MD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-MD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(mdM(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(mdM(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                % EVE2
                for j=1:length(EVE2labelNames)
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-MD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-MD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(md1(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(md1(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-MD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-MD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(md2(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(md2(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-MD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-MD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(mdM(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(mdM(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                % EVE3
                for j=1:length(EVE3labelNames)
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-MD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-MD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(md1(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(md1(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-MD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-MD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(md2(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(md2(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-MD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-MD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(mdM(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(mdM(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                % SLANT BC
                for j=1:length(BClabelNames)
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-MD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-MD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(md1(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(md1(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-MD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-MD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(md2(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(md2(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-MD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-MD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(mdM(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(mdM(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                %%%%%%%%%%%%%%% AD
                for j=1:length(EVE1labelNames)
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-AD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-AD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(ad1(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(ad1(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-AD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-AD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(ad2(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(ad2(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-AD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-AD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(adM(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(adM(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                for j=1:length(EVE2labelNames)
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-AD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-AD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(ad1(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(ad1(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-AD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-AD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(ad2(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(ad2(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-AD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-AD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(adM(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(adM(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                for j=1:length(EVE3labelNames)
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-AD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-AD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(ad1(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(ad1(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-AD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-AD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(ad2(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(ad2(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-AD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-AD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(adM(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(adM(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                for j=1:length(BClabelNames)
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-AD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-AD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(ad1(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(ad1(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-AD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-AD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(ad2(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(ad2(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-AD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-AD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(adM(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(adM(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                
                %%%%%%%%%%%%%%% RD
                for j=1:length(EVE1labelNames)
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-RD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI1-RD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(rd1(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(rd1(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-RD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTI2-RD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(rd2(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(rd2(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-RD-mean'];
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'DTIM-RD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(rdM(eve1(:)==EVE1labelID(j)));
                        ColValues{end+1} = std(rdM(eve1(:)==EVE1labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                for j=1:length(EVE2labelNames)
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-RD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI1-RD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(rd1(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(rd1(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-RD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTI2-RD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(rd2(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(rd2(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-RD-mean'];
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'DTIM-RD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(rdM(eve2(:)==EVE2labelID(j)));
                        ColValues{end+1} = std(rdM(eve2(:)==EVE2labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                for j=1:length(EVE3labelNames)
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-RD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI1-RD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(rd1(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(rd1(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-RD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTI2-RD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(rd2(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(rd2(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-RD-mean'];
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'DTIM-RD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(rdM(eve3(:)==EVE3labelID(j)));
                        ColValues{end+1} = std(rdM(eve3(:)==EVE3labelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                for j=1:length(BClabelNames)
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-RD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI1-RD-std'];
                    if DTI1 > 0
                        ColValues{end+1} = mean(rd1(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(rd1(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-RD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTI2-RD-std'];
                    if DTI2 > 0
                        ColValues{end+1} = mean(rd2(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(rd2(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                    
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-RD-mean'];
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'DTIM-RD-std'];
                    if DTIM > 0
                        ColValues{end+1} = mean(rdM(bc(:)==BClabelID(j)));
                        ColValues{end+1} = std(rdM(bc(:)==BClabelID(j)));
                    else
                        ColValues{end+1} = NaN;
                        ColValues{end+1} = NaN;
                    end
                end
                
                %%%%%%%%%%%%%%% ROI Volumes
                for j=1:length(BClabelNames)
                    ColHeader{end+1} = ['BrainColor-' BClabelNames{j} '-' 'Volume'];
                    ColValues{end+1} = sum(bc(:)==BClabelID(j))*prod(bcinfo.PixelDimensions(1:3));
                end
                
                for j=1:length(EVE1labelNames)
                    ColHeader{end+1} = ['EveType1-' EVE1labelNames{j} '-' 'Volume'];
                    ColValues{end+1} = sum(eve1(:)==EVE1labelID(j))*prod(eve1info.PixelDimensions(1:3));
                end
                for j=1:length(EVE2labelNames)
                    ColHeader{end+1} = ['EveType2-' EVE2labelNames{j} '-' 'Volume'];
                    ColValues{end+1} = sum(eve2(:)==EVE2labelID(j))*prod(eve2info.PixelDimensions(1:3));
                end
                for j=1:length(EVE3labelNames)
                    ColHeader{end+1} = ['EveType3-' EVE3labelNames{j} '-' 'Volume'];
                    ColValues{end+1} = sum(eve3(:)==EVE3labelID(j))*prod(eve3info.PixelDimensions(1:3));
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