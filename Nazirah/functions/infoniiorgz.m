function info = infoniiorgz(filename)


try
    if strcmpi(filename(end-1:end),'gz')
        if ~isfile(filename)
            filename = filename(1:end-3); % remove .gz to load
        end
    elseif strcmpi(filename(end-2:end),'nii')
        if ~isfile(filename)
            filename = [filename '.gz']; % add .gz to load
        end
    end
    
    %fprintf('Loading filename: %s\n', filename)
    info = niftiinfo(filename);
catch err
    info = [];
    fprintf('ERROR: %s\n',err.message)
end
