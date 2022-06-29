function [bounding_box,mask_volume] = find_bounding_box(label_volume,labels,pad_distance,dilate_distance)
%FIND_BOUNDING_BOX - Finds the bounding box givens a set of labels
%
% Syntax:  [bounding_box,mask_volume] = find_bounding_box(input1,input2,input3)
%
% Inputs:
%    label_volume    - Volume to find the mask in
%    labels          - labels to include in mask
%    pad_distance    - distance to pad bounding box
%    dilate_distance - distance to dilate mask
%
% Outputs:
%    bounding_box    - 3x2 matrix containing the bouning box of the labels
%    mask_volume     - mask of where the labels occur
%
% See also: extract_bounding_box
%
% Author:  plassaaj
% Date:    13-Feb-2015 12:30:13
% Version: 1.0
% Changelog:
%
% 13-Feb-2015 12:30:13 - initial creation
%
%------------- BEGIN CODE --------------
    mask_volume = zeros(size(label_volume));
    for i=1:numel(labels)
        l = labels(i);
        mask_volume = mask_volume + (label_volume == l);
    end
    
    bounding_box = zeros(3,2);
    % x bounding
    dim = 1;
    for i=1:size(mask_volume,dim)
        slice = mask_volume(i,:,:);
        if any(slice(:))
            bounding_box(dim,1) = max([i-pad_distance 1]);
            break
        end
    end
    for i=size(mask_volume,dim):-1:1
        slice = mask_volume(i,:,:);
        if any(slice(:))
            bounding_box(dim,2) = min([i+pad_distance size(mask_volume,dim)]);
            break
        end
    end
    % y bounding
    dim = 2;
    for i=1:size(mask_volume,dim)
        slice = mask_volume(:,i,:);
        if any(slice(:))
            bounding_box(dim,1) = max([i-pad_distance 1]);
            break
        end
    end
    for i=size(mask_volume,dim):-1:1
        slice = mask_volume(:,i,:);
        if any(slice(:))
            bounding_box(dim,2) = min([i+pad_distance size(mask_volume,dim)]);
            break
        end
    end
    % z bounding
    dim = 3;
    for i=1:size(mask_volume,dim)
        slice = mask_volume(:,:,i);
        if any(slice(:))
            bounding_box(dim,1) = max([i-pad_distance 1]);
            break
        end
    end
    for i=size(mask_volume,dim):-1:1
        slice = mask_volume(:,:,i);
        if any(slice(:))
            bounding_box(dim,2) = min([i+pad_distance size(mask_volume,dim)]);
            break
        end
    end
    mask_volume = imdilate(mask_volume,ones(dilate_distance,dilate_distance,dilate_distance));
end

