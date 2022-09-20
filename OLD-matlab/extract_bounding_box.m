function [new_image] = extract_bounding_box(img,bounding_box)
%EXTRACT_BOUNDING_BOX - Reduces an image to it's bounding box
%
% Syntax:  [new_image] = extract_bounding_box(image,bounding_box)
%
% Inputs:
%    img          - 3D input image volume
%    bounding_box - 3x2 bounding box of the input image
%
% Outputs:
%    new_image    - img reduced to bounding box
%
% See also: find_bounding_box
%
% Author:  plassaaj
% Date:    13-Feb-2015 12:43:33
% Version: 1.0
% Changelog:
%
% 13-Feb-2015 12:43:33 - initial creation
%
%------------- BEGIN CODE --------------

new_image = img(bounding_box(1,1):bounding_box(1,2),...
    bounding_box(2,1):bounding_box(2,2),...
    bounding_box(3,1):bounding_box(3,2));
end