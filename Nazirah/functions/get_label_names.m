function [EVElabelNames,EVElabelID,BClabelNames,BClabelID,DTISeglabelNames,DTISeglabelID] = get_label_names(EVE_path,BC_path,DTI_path)


EVElabelNames = readcell(EVE_path,"delimiter",'\n');
for i=1:length(EVElabelNames)
    EVElabelID(i) = sscanf(EVElabelNames{i},'%d');
end

BClabelNames = readcell(BC_path,"delimiter",'\n');
BClabelNames = BClabelNames(2:end); %remove header
BClabelNames2 = split(BClabelNames,',');
BClabelNames = join(BClabelNames2(:,[2 1]),' ');
for i=1:length(BClabelNames2)
    BClabelID(i) = str2double(BClabelNames2{i,2});
end

DTISeglabelNames = readcell(DTI_path,"delimiter",'\n');
for i=1:length(DTISeglabelNames)
    DTISeglabelID(i) = sscanf(DTISeglabelNames{i},'%d',1);
end