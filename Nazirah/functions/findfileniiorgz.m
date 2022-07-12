function filename = findfileniiorgz(filedir,filestring)

% remove filesep in filedir first
if strcmp(filedir(end),'/') || strcmp(filedir(end),'\')
    filedir(end) = [];
end

% find file - should be either nii or nii.gz
filedata = dir([filedir filesep '*' filestring '*']);

if length(filedata) > 1
    fprintf('Warning.. there are multiple files. Choosing the first\n')
end

% get the fullpath of the file
filename = [filedata(1).folder filesep filedata(1).name];