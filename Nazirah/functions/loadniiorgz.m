function data = loadniiorgz(filename)

try
    if strcmpi(filename(end-1:end),'gz')
        if isfile(filename)
            data = niftiread(filename);
        else
            data = niftiread(filename(1:end-2));
        end
    elseif strcmpi(filename(end-2:end),'nii')
        if isfile(filename)
            data = niftiread(filename);
        else
            data = niftiread([filename '.gz']);
        end
    end
catch err
    data = [];
    fprintf('ERROR: %s\n',err.message)
end