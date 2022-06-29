function raw2nii(in_raw, out_nii, datatype, dims, res, ort)

% raw2nii('U:\CV\Users\SRI\VUDBS_MATLAB_PROGRAMS\RAW2NIFTI\Test_Raw_to_NIFTI.img', 'U:\CV\Users\SRI\VUDBS_MATLAB_PROGRAMS\RAW2NIFTI\Test_Raw_to_NIFTI.nii', 'int16', [256 256 190], [1 1 1], {'AP', 'SI', 'RL'});
% raw2nii('U:\CV\Users\SRI\VUDBS_MATLAB_PROGRAMS\RAW2NIFTI\Test_Raw_to_NIFTI.img', 'U:\CV\Users\SRI\VUDBS_MATLAB_PROGRAMS\RAW2NIFTI\Test_Raw_to_NIFTI.nii', 'int16', [256 256 190], [1 1 1], {2, -3, 1});

% raw2nii - function that converts raw image format to a nifti file
%
% Input: in_raw - the input raw file
%        out_nii - the output nii file
%        datatype - the datatype of the image (string)
%                 - e.g., 'single', 'double', 'uint16', 'int8', 'int16'
%        dims - the dimensions of the input image (1 x 3 vector)
%             - e.g., [256 256 312]
%        res - the resolution of the input image (1 x 3 vector)
%            - e.g., [0.5 0.5 3]
%            - assumed to be in millimeters
%            - the order of "res" should correspond to the order of "dims"
%        ort - the orientation of the input image (1 x 3 cell array)
%            - e.g., {'RL', 'SI', 'PA'} would indicate that:
%              -> the first dimension is oriented Right to Left
%              -> the second dimension is oriented Superior to Inferior
%              -> the third dimension is oriented Posterior to Anterior
%            - the only valid values for "ort" are:
%              -> 'RL' -- Right to Left
%              -> 'LR' -- Left to Right
%              -> 'SI' -- Superior to Inferior
%              -> 'IS' -- Inferior to Superior
%              -> 'PA' -- Posterior to Anterior
%              -> 'AP' -- Anterior to Posterior
%
% Output: None, however, the nii file is saved as "out_nii"
%
% Written by: Andrew Asman (5/10/2013), Modified by Yuann Liu (6/5/2013), Modified by Srivatsan Pallavaram (1/29/2015)
% Example command call:  raw2nii('T1.im', 'T1.nii', 'uint16', [256 256 170], [1 1 1], {'AP', 'SI', 'RL'})


addpath('U:\CV\Users\SRI\VUDBS_MATLAB_PROGRAMS\niitools\');
% read in the raw file ()
fid = fopen(in_raw, 'r');
img = fread(fid, prod(dims), datatype);
fclose(fid);

% reshape to the proper dimensions
img = reshape(img, dims);

% fix the order of the dimensions
order = zeros([1 3]);
 % for k = 1:3
    %     if (strcmp(ort{k}, 'LR') || strcmp(ort{k}, 'RL'))
    %         order(1) = k;
    %     elseif (strcmp(ort{k}, 'PA') || strcmp(ort{k}, 'AP'))
    %         order(2) = k;
    %     elseif (strcmp(ort{k}, 'SI') || strcmp(ort{k}, 'IS'))
    %         order(3) = k;
    %     else
    %         error(sprintf('Unrecognized orientation %s', ort{k}));
    %     end
    % end
    
    for k = 1:3
        if (strcmp(ort{k}, 'LR') || strcmp(ort{k}, 'RL') || (ort{k} == 1) ||  (ort{k} == -1))
            order(1) = k;
        elseif (strcmp(ort{k}, 'PA') || strcmp(ort{k}, 'AP') || (ort{k} == 2) ||  (ort{k} == -2))
            order(2) = k;
        elseif (strcmp(ort{k}, 'SI') || strcmp(ort{k}, 'IS') || (ort{k} == 3) ||  (ort{k} == -3))
            order(3) = k;
        else
            error(sprintf('Unrecognized orientation %s', ort{k}));
        end
    end
    
    dims = dims(order);
    res = vox(order);
    ort = ort(order);
    img = permute(img, order);
    
    % flip any dimensions that we need to
    if (strcmp(ort{1}, 'RL') || (ort{1} == 1))
        img = flipdim(img, 1);
    end
    if (strcmp(ort{2}, 'AP') || (ort{2} == 2))
        img = flipdim(img, 2);
    end
    if (strcmp(ort{3}, 'SI') || (ort{3} == -3))
        img = flipdim(img, 3);
    end
    
% set the datatype
dt = 16;
switch datatype
    case 'uint8'
        dt = 2;
    case 'int16'
        dt = 4;
    case 'int32'
        dt = 8;
    case 'single'
        dt = 16;
    case 'double'
        dt = 64;
    case 'int8'
        dt = 256;
    case 'uint16'
        dt = 512;
    case 'uint32'
        dt = 768;
    otherwise
        dt = 16;
end

% create the nii file
nii = make_nii(img, res, [0 0 0], dt);
nii.hdr.dime.xyzt_units = 10;
nii.hdr.dime
nii.hdr.hist

% save the nii file
save_nii(nii, out_nii);

