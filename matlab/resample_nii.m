function resample_nii(in_fname, out_fname, out_res, interp_type);
% RESAMPLE_NII - resamples an input nifti file
%
% resample_nii(in_fname, out_fname, out_res, interp_type);
%
% Input: in_fname - the input nifti file
%        out_fname - the resampled nifti file
%        out_res - the output voxel resolution (e.g., [2 2 2]);
%        interp_type - the interpolation type to use (e.g., 'linear')
%                    - see "help interpn" for the available options
%
% Output: (None)
%

a = [];
typ = 0;
if length(in_fname) >= 7 && strcmp(in_fname(end-6:end), '.nii.gz')
    a = load_untouch_nii_gz(in_fname);
    typ = 0;
elseif length(in_fname) >= 4 && strcmp(in_fname(end-3:end), '.nii');
    a = load_untouch_nii(in_fname);
    typ = 1;
else
    error('invalid input file -- must .nii.gz or .nii');
end

if length(out_res) ~= 3
    error('out_res should be a vector of the form [Rx Ry Rz]');
end

% get the info from the nifti file
res = a.hdr.dime.pixdim(2:4);
dims = a.hdr.dime.dim(2:4);
rates = out_res ./ res;

% do the resampling
if dims(3) > 1

    % set up the grid
    [x1 y1 z1] = ndgrid(1:dims(1), ...
                        1:dims(2), ...
                        1:dims(3));
    [x2 y2 z2] = ndgrid(1:rates(1):dims(1), ...
                        1:rates(2):dims(2), ...
                        1:rates(3):dims(3));

    % modify the nifti file
    a.img = interpn(x1, y1, z1, double(a.img), x2, y2, z2, interp_type);
    a.hdr.dime.dim(2:4) = size(a.img);
    a.hdr.dime.pixdim(2:4) = out_res(1:3);
else

    % set up the 2D grid
    [x1 y1] = ndgrid(1:dims(1), ...
                     1:dims(2));
    [x2 y2] = ndgrid(1:rates(1):dims(1), ...
                     1:rates(2):dims(2));

    % modify the nifti file
    a.img = interpn(x1, y1, double(a.img), x2, y2, interp_type);
    a.hdr.dime.dim(2:3) = size(a.img);
    a.hdr.dime.pixdim(2:3) = out_res(1:2);
end

% save the result
if (typ == 0)
    save_untouch_nii_gz(a, out_fname);
else
    save_untouch_nii(a, out_fname);
end

