function results = compareFiles(adpath,rdpath)

%find folder
adpath(strfind(adpath,'ad2.nii'):end) = [];
rdpath(strfind(rdpath,'rd2.nii'):end) = [];

try
    ad1 = niftiread([adpath 'ad1.nii.gz']);
catch
    ad1 = niftiread([adpath 'ad1.nii']);
end

try
    rd1 = niftiread([rdpath 'rd1.nii.gz']);
catch
    rd1 = niftiread([rdpath 'rd1.nii']);
end

try
    ad2 = niftiread([adpath 'ad2.nii.gz']);
catch
    ad2 = niftiread([adpath 'ad2.nii']);
end

try
    rd2 = niftiread([rdpath 'rd2.nii.gz']);
catch
    rd2 = niftiread([rdpath 'rd2.nii']);
end

%comparepath = strsplit(adpath,'/');

%results.session = comparepath{end-4};
results.ad = isequal(ad1,ad2);
results.rd = isequal(rd1,rd2);

end

