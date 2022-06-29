function nii = loadniiorniigz( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(length(dir([filename '.gz']))>0)
    
    nii = load_nii_gz([filename '.gz']);
else
    if(length(dir(filename))>0)
        if(and(filename(end-1)=='g',filename(end)=='z'))
            nii=load_untouch_nii_gz(filename);
        else
            nii = load_untouch_nii(filename);
        end
    else
        disp('cannot find file');
        throw(['cannot file file: ' filename]);
    end
    
end
