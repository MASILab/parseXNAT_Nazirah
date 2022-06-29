function [img] = masi_load_image(filename)
%Generic image loading function which loads an image of any type
%The following formats are supported:
%   Nii
%   Nii.gz
%   par/rec
%   dicom (.dcm,.ima)
%   ALSO Supported: .mat and non-medical image types: png,jpg,tiff,etc...
%
%INPUTS:
%   filename - the full path to the image to be loaded
%
%OUTPUTS:
%   img - struct - image which was loaded. data will always be in
%   img.image_data, which may be a cell array if there are multiple volumes
%   all other fields are any other header information which
%   could be loaded
%
%Author: Rob Harrigan
%
%Changelog:
%   Rob Harrigan - 11/25/2014 - initial commit
%
%Functionality to be added:
%   rec.gz
%   .img/.hdr files
%   load dicominfo?

[path,name,ext] = fileparts(filename);

if ~exist(filename,'file')
    error('File does not exist');
end

switch lower(ext)
    case '.nii'
        img = load_nii(filename); %untested
        img.image_data = img.img;
        img = rmfield(img,'img');
    case '.gz'
        split_name = strsplit(name,'.');
        ext2 = split_name{end};
        switch lower(ext2)
            case 'nii'
                img = load_nii_gz(filename);
                img.image_data = img.img;
                img = rmfield(img,'img');
            otherwise
                error('Unsupported file type');
        end
    case '.dcm'
        vol = dicomread(filename);
        img.image_data = vol;
    case '.ima'
        vol = dicomread(filename); %WARNING untested - no files available
        img.image_data = vol;
    case '.par'
        par = loadPAR(filename);
        name = strrep(name,'PAR','REC');
        name = strrep(name,'par','rec');
        rec_file = fullfile(path,[name,'.rec']);
        cap_rec_file = fullfile(path,[name,'.REC']);
        rec_gz_file = fullfile(path,[name,'.rec.gz']);
        if exist(rec_file,'file')
            rec = loadREC(rec_file,par,[],false);
        elseif exist(cap_rec_file,'file')
            rec = loadREC(cap_rec_file,par,[],false);
        elseif exist(rec_gz_file,'file')
            error('Unsupported file type');
        else
            error('No corresponding rec file with provided par file.');
        end
        img = par;
        img.image_data = rec;
    case '.rec'
        name = strrep(name,'REC','PAR');
        name = strrep(name,'rec','par');
        par_file = fullfile(path,[name,'.par']);
        cap_par_file = fullfile(path,[name,'.PAR']);
        if exist(par_file,'file')
            par = loadPAR(par_file);
        elseif exist(cap_par_file,'file')
            par = loadPAR(cap_par_file);
        else
            error('No corresponding par file with provided rec file. exiting');
        end
        rec = loadREC(filename,par,[],false);
        img = par;
        img.image_data = rec;
    case '.hdr'
        error('Unsupported file type');
    case '.img'
        error('Unsupported file type');
    otherwise 
        try
            im = imread(filename);
            img.image_data = im;
        catch
            warning('Unsupported file type, please read help, exiting');
            exit
        end
end
    


end