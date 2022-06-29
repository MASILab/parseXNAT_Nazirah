function [data] = apply_smooth_boundary_model(truth, shiftiness, span);
% APPLY_BOUNDARY_MODEL - Applies the gaussian-based shifting model to a given
%   set of truth labels.
%
% [data] = apply_boundary_model(truth, shiftiness);
%
% Input: truth - the true "segmentation"
%        shiftiness - the amount of error exhibited at each boundary
% Output: data - the observation of the true segmentation
%

dims = size(truth);
if length(dims) == 2
    dims = [dims 1];
end
data = truth;

% only perform the shifting on the x and y components

xshift_list_t = zeros(500000, 5);
yshift_list_t = zeros(500000, 5);

xs = 1:dims(1)-1;
xn = 2:dims(1);

ys = 1:dims(2)-1;
yn = 2:dims(2);

% get the xshift_list and yshift_list
cx = 0;
cy = 0;
for z = 1:dims(3)
    for y = 1:dims(2);
        for x = 1:(dims(1)-1)

            if (x ~= dims(1))
                if (truth(x, y, z) ~= truth(x+1, y, z))
                    cx = cx + 1;
                    vals = [truth(x, y, z), truth(x+1, y, z)];
                    if x < dims(1)/2
                        xshift_list_t(cx, :) = [x+1, y, z, vals(1), vals(2)];
                    else
                        xshift_list_t(cx, :) = [x, y, z, vals(1), vals(2)];
                    end
                end
            end

            if (y ~= dims(2))
                if (truth(x, y, z) ~= truth(x, y+1, z))
                    cy = cy + 1;
                    vals = [truth(x, y, z), truth(x, y+1, z)];
                    if y < dims(2)/2
                        yshift_list_t(cy, :) = [x, y+1, z, vals(1), vals(2)];
                    else
                        yshift_list_t(cy, :) = [x, y, z, vals(1), vals(2)];
                    end
                end
            end
        end
    end
end
xshift_list = xshift_list_t(1:cx, :);
yshift_list = yshift_list_t(1:cy, :);

% randomly generate the x and y shifts
xshifts = randn([size(xshift_list, 1) 1]) * shiftiness;
yshifts = randn([size(yshift_list, 1) 1]) * shiftiness;

% introduce the spatial correlations (but maintain the same magnitude)
xshifts = round(smooth(xshifts, span));
yshifts = round(smooth(yshifts, span));

% apply the x-shifts
for i = 1:size(xshift_list, 1)
    x = xshift_list(i, 1);
    y = xshift_list(i, 2);
    z = xshift_list(i, 3);
    v1 = xshift_list(i, 4);
    v2 = xshift_list(i, 5);
    shift = xshifts(i);

    if shift < 0
        data(max([x+shift 1]):x, y, z) = v2;
    elseif shift > 0
        data(x:min([x+shift dims(1)]), y, z) = v1;
    end
end

% apply the y-shifts
for i = 1:size(yshift_list, 1)
    x = yshift_list(i, 1);
    y = yshift_list(i, 2);
    z = yshift_list(i, 3);
    v1 = yshift_list(i, 4);
    v2 = yshift_list(i, 5);
    shift = yshifts(i);

    if shift < 0
        data(x, max([y+shift 1]):y, z) = v2;
    elseif shift > 0
        data(x, y:min([y+shift dims(2)]), z) = v1;
    end
end

