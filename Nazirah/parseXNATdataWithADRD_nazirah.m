

% Setup the root directory for all subjects
D = '/fs4/masi/landmaba/BLSAdti/BLSA/'

addpath('/fs4/masi/landmaba/BLSAdti/BLSA')
addpath(genpath('/fs4/masi/landmaba/BLSAdti/matlab'))
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

for jSubj = 1:length(SUBJS) %929 % [73 128 450 528 544 569 586 621 720 834 929]%834% 1:length(SUBJS)
    
    % Get a list of Sessions; ignore . and ..
    %SESSIONS = dir([D SUBJS(jSubj).name]); SESSIONS =SESSIONS(3:end);
    SESSIONS = dir([D SUBJS(jSubj).name]);
    SESSIONS(startsWith({SESSIONS.name},'.')) = [];
    
    for jSession=1:length(SESSIONS)     
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

                
                %% Deal with the first DTI session "DTI(1)"
                dti1
                % if ver3 == 1, don't run verifyADRD
                
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
                label1name = [DS Stamper(1).name filesep 'WM_LABELS' '/Rectified_EVE_Labels.nii.gz'];
                eveName = [label1name '.subjLabels.nii.gz' ];
                if(length(dir(eveName))<1)
                    % reg_transform -ref ref_img -invAff transform_mat transform_mat.inv
                    % link: http://cmictig.cs.ucl.ac.uk/wiki/index.php/Reg_transform
                    system(['reg_transform -ref ' fa1name ' -invAff ' xfm1name ' ' xfm1name '.inv']);
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
                
                faM = loaduntouchniiorniigz(faMname);
                faM2 = niftiread(faMname);
                
                
            end
        catch err
            disp('Oh no!')
            fp=fopen([D 'FatalErrorsWithADRDv8.txt'],'at');
            fprintf(fp,'%s - %s\n',SESSIONS(jSession).name,err.message);
            fclose(fp);
        end
    end
end