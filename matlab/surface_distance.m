function [msd mesd hd] = surface_distance(truth, estimate, res_dims)
% SURFACE_DISTANCE - returns the surface distance metrics between a true
%                    segmentation and an estimate
%
% [msd mesd hd] = surface_distance(truth, estimate)
%
% Input: truth - the true segmentation
%        estimate - the estimated segmentation
%        res_dims - the resolution of each of the dimensions
%                   should be an array of the form [Rx Ry Rz]
%
% Output: msd - the mean surface distance
%         mesd - the medisn surface distance (not effected by outliers like mean)
%         hd - the hausdorff distance (max distance)


% convert the truth and estimate to isosurfaces
iso1 = isosurface(truth == 1, 0);
iso2 = isosurface(estimate == 1, 0);
vert1 = iso1.vertices;
vert2 = iso2.vertices;
clear iso1 iso2

% get the number of elements in each of the isosurfaces
num1 = size(vert1, 1);
num2 = size(vert2, 1);

% store the distance array
dist = zeros([num1, 1]);

% allocate some temporary matrices to make it faster
rd = repmat(res_dims, [num2 1]);
om = ones(num2, 1);

% calculate the distance from the estimate to the truth for each element
for i = 1:num1
    dist(i) = sqrt(min(sum((rd .* (vert2 - om*vert1(i, :))).^2, 2)));
end

% get the mean surface distance (msd) and the hausdorff distance (hd)
msd = mean(dist);
mesd = median(dist);
hd = max(dist);

