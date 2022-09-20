

% Setup the root directory for all subjects
%D = '/fs4/masi/landmaba/BLSAdti/BLSA/'

%addpath('/fs4/masi/landmaba/BLSAdti/BLSA')
%addpath(genpath('/fs4/masi/landmaba/BLSAdti/matlab'))

%D = '/Users/nana/Documents/MATLAB/BLSA/original_files/Nazirah/';

D = '/home/local/VANDERBILT/mohdkhn/Documents/BLSA_original/data';
addpath('/home/local/VANDERBILT/mohdkhn/Documents/BLSA/Nazirah/functions');
% this m-file should be in MASI-42 /home/local/VANDERBILT/mohdkhn/Documents/BLSA/Nazirah

%%
% Load the label names
EVE_path = [D 'EVE_Labels.csv'];
BC_path = [D 'andrew_multiatlas_labels.csv'];
DTI_path = [D 'DTIseg.txt'];

[EVElabelNames,EVElabelID,BClabelNames,BClabelID,DTISeglabelNames,DTISeglabelID] = get_label_names(EVE_path,BC_path,DTI_path);

%% Get a list of Subjects; ignore . and ..
%SUBJS = dir([D filesep 'BLSA*']); SUBJS=SUBJS(1:end);
SUBJS = dir([D 'BLSA*']);


%% This is used to ensure that the stats file is replaced each time.
doMakePNG = 0;
doSkipPNGDone = 0;
doMakeReport = 1;
results = [];

for jSubj = 1:length(SUBJS) %929 % [73 128 450 528 544 569 586 621 720 834 929]%834% 1:length(SUBJS)
    
    % Get a list of Sessions; ignore . and ..
    %SESSIONS = dir([D SUBJS(jSubj).name]); SESSIONS =SESSIONS(3:end);
    SESSIONS = dir([D SUBJS(jSubj).name]);
    SESSIONS(startsWith({SESSIONS.name},'.')) = [];
    
    for jSession=1%:length(SESSIONS)     
        DS = [D SUBJS(jSubj).name filesep SESSIONS(jSession).name filesep];
        Stamper = dir([DS  '*Stamper*']);
        
        ver3 = 0;
        dtiQA = dir([DS '*dtiQA_v2_1']);
        if(length(dtiQA)~=2)
            dtiQA = dir([DS '*dtiQA_v2']);
            disp('Got dtiQA_v2')
        end
        
        if(length(dtiQA)~=2)
            dtiQA = dir([DS '*dtiQA_v2*']);
            disp('Got *dtiQA_v2*')
        end
        
        if((length(dtiQA)~=2)+strcmp(SESSIONS(jSession).name,'BLSA_4889_02-0_10')+strcmp(SESSIONS(jSession).name,'BLSA_4806_05-0_10'))
            ver3 = 1;
            dtiQA = dir([DS '*dtiQA_v3*']);
            disp('Got *dtiQA_v3*')
        end
        dtiQAMulti = dir([DS '*dtiQA_Multi*']);
        MultiAtlas = dir([DS '*Multi_Atlas*']);
        
        % Sometimes Multi-atlas will find another T1. If
        % more than one are found, choose the MPRAGE one.
        if(length(MultiAtlas)>1)
            MultiAtlas = dir([DS '*MPRAGE*Multi_Atlas*']);
            disp('More than 1 multiatlas')
        end
        MPRAGE = dir([DS 'MPRAGE*']);
        
        disp([SUBJS(jSubj).name ' ' SESSIONS(jSession).name])
        disp([length(Stamper) length(dtiQA) length(MultiAtlas)])
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
                    reportFileName = [D 'statsWithADRDVol' filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
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
                disp('FOUND valid dataset')
                
                % Create CSV report file.
                reportFileName = [D 'statsWithADRDVol' filesep SESSIONS(jSession).name '-AllStatsWithADRDVol.csv'];
                
                % If CSV file exist, don't continue below. Move to the next Session
                if(exist(reportFileName,'file'))
                    disp(['Already done: ' reportFileName]);
                    continue;
                end
                disp('touch')
                system(['echo `date` > ' reportFileName])
                
                % If PNG exists, don't continue below and go to next session. [Whyy?]
                if(doSkipPNGDone)
                    if(exist([D 'pngs' filesep SESSIONS(jSession).name '.png'],'file'))
                        disp('It''s done');
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
                
                %% REMOVE AFTER TEST
                multiDTIold %save files as ad2,rd2
                
                results(end+1).session = [dtiQAMulti(1).name];
                results(end).results_compare = compareFiles(adMname,rdMname);
                
                %% Deal with the first DTI session "DTI(1)"
                dti1
                % if ver3 == 1, don't run verifyADRD
                %% REMOVE AFTER TEST
                dti1old
                
                results(end+1).session = [dtiQA(1).name];
                results(end).results_compare = compareFiles(ad1name,rd1name);
                
                % DELETE all test files
                adrdpath = [DS dtiQA(1).name filesep 'AD' filesep];
                delete([adrdpath 'ad1.nii'],[adrdpath 'rd1.nii'],[adrdpath 'ad2.nii.gz'],[adrdpath 'rd2.nii.gz']);
            end
            
        catch err
            disp('Oh no!')
            fp=fopen([D 'FatalErrorsWithADRDv8.txt'],'at');
            fprintf(fp,'%s - %s\n',SESSIONS(jSession).name,err.message);
            fclose(fp);
        end
    end
end
%% REMOVE AFTER TEST - save test results
save(sprintf('compareADRD_nii_%s',date),results)