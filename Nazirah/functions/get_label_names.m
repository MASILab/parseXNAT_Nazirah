function [EVElabelNames,EVElabelID,BClabelNames,BClabelID,DTISeglabelNames,DTISeglabelID] = get_label_names(EVE_path,BC_path,DTI_path)


EVElabelNames = readcell(EVE_path,"delimiter",'\n');
for i=1:length(EVElabelNames)
    EVElabelID(i) = sscanf(EVElabelNames{i},'%d');
end

BClabelNames = readcell(BC_path,"delimiter",'\n');
for i=1:length(BClabelNames)
    BClabelID(i) = sscanf(BClabelNames{i},'%d ',1);
end

DTISeglabelNames = readcell(DTI_path,"delimiter",'\n');
for i=1:length(DTISeglabelNames)
    DTISeglabelID(i) = sscanf(DTISeglabelNames{i},'%d',1);
end