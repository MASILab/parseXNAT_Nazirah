function colored_boxplot(data3d, varargin)
% colored_boxplot - generates a colored boxplot for comparing accuracy of
%                   multiple algorithms across multiple categories
%                   (e.g., labels)
%
% Two forms:
% 1) colored_boxplot(data3d)
% 2) colored_boxplot(data3d, opts)
%
% Input: data3d - should be an X x Y x Z matrix, where:
%               - X is the number of samples
%               - Y is the number of labels
%               - Z is the number of algorithms
%        opts - the options struct
%
% Output: None
%
% *** The options struct ***
%
% opts.ylimits - the limits for the y axis -- e.g., [0 1]
% opts.colors - the colors for each algorithm (default: hsv)
% opts.label_names - the names for each label (cell array) (default: none)
% opts.names - the names for each algorithm (cell array) (default: none)

if length(varargin) == 0
    opts = struct;
elseif length(varargin) == 1
    opts = varargin{1};
else
    error('Too many input arguments');
end

if length(size(data3d)) == 2
    data3d = permute(data3d, [1 3 2]);
end

% derived settings
num_labs = size(data3d, 2);
num_algs = size(data3d, 3);

% convert the 3d data into a 2D format
data2d = zeros(size(data3d, 1), size(data3d, 2)*size(data3d, 3));
cc = 1;
for l = 1:num_labs
    for e = 1:num_algs
        data2d(:, cc) = data3d(:, l, e);
        cc = cc + 1;
    end
end

if ~isfield(opts, 'ylimits')
    opts.ylimits(1) = min(data2d(:)) - std(data2d(:)) / 2;
    opts.ylimits(2) = max(data2d(:)) + std(data2d(:)) / 2;
end

% create the initial boxplot
h = boxplot(data2d);

% add the colors
if ~isfield(opts, 'colors')
    opts.colors = hsv(num_algs);
end
for k = 1:(num_labs*num_algs)
    set(h(:, k), 'Color', opts.colors(1 + mod(k - 1, num_algs), :), ...
                 'LineWidth', 2);
end

% add the separating lines
for l = 1:(num_labs-1)
    hold on;
    xx = l*num_algs + 0.5;
    plot([xx xx], opts.ylimits, 'k', 'LineWidth', 2);
    hold off;
end

if isfield(opts, 'label_names')
    % add the label names
    for l = 1:num_labs
        xx = l*num_algs - ((num_algs-1)/2);
        text(xx, opts.ylimits(2)+0.01, opts.label_names{l}, ...
             'HorizontalAlignment', 'center');
    end
end

% set the ylimits
ylim(opts.ylimits);

if isfield(opts, 'names')
    % add the legend
    legend(h(3, 1:num_algs), opts.names{:}, 'Location', 'Best');
end

grid minor;
grid on;
set(gca, 'XTickLabel', {''});


