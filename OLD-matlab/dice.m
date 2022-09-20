function [meandic dic uniquelabels] = dice(truth, estimate, varargin)
% DICE - calculate the dice similarity coefficient
%
% [meandic dic uniquelabels] = dice(truth, estimate)
%
% Assumes: truth, estimate are the same dimension, orientation, etc...
%          only interested in the labels that are unique to the truth
% Input: truth - truth volume/vector/matrix
%        estimate - estimate volume/vector/matrix
%
% Output: meandic - mean DSC excluding the first label (presumably background)
%         dic - DSC value for each label
%         uniquelabels - the label numbers corresponding to 'dic'

if length(varargin) == 1
    list = varargin{1};
    for i = 1:length(list)
        truth(truth == list(i)) = 0;
        estimate(estimate == list(i)) = 0;
    end
end

% find the unique values in the truth
uniquelabels = unique(truth);

% calculate the dice similarity coefficient (DSC) for each label
dic = zeros([length(uniquelabels) 1]);

% iterate over each unique value
for i = 1:length(uniquelabels)

    % get the curent label of interest
    v = uniquelabels(i);

    % set the binary images for this label
    tl = truth == v;
    el = estimate == v;

    % get the number non-zeros (nnz) for each image and their intersect
    nnz_tl = nnz(tl);
    nnz_el = nnz(el);
    nnz_and = nnz(tl&el);

    % set the dice value
    if (nnz_tl == 0 & nnz_el == 0)
        dic(i) = 1;
    else
        dic(i) = 2*nnz_and / (nnz_tl + nnz_el);
    end
end

meandic = mean(dic(2:end)); %ignores background

